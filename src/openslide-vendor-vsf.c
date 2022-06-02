/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2007-2014 Carnegie Mellon University
 *  Copyright (c) 2011 Google, Inc.
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, version 2.1.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with OpenSlide. If not, see
 *  <http://www.gnu.org/licenses/>.
 *
 */

/*
 * VSF image format
 *
 */

#include <config.h>

#include <string.h>
#include <sys/types.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <glib.h>

#include "openslide-private.h"
#include "openslide-decode-gdkpixbuf.h"
#include "openslide-decode-jpeg.h"
#include "openslide-decode-jp2k.h"
#include "openslide-decode-png.h"
#include "openslide-decode-xml.h"

//#include "openslide-vendor-vsf.h"

#define INDEX_FILE_EXTENSION ".vsf"
#define IMAGE_FILE_EXTENSION ".img"
#define VMSCOPE_PRODUCT_HEADER "VSF%c.%c VMscope GmbH (Germany)"
#define OPENSLIDE_PROPERTY_VSF_FILENAME "vsf.filename"

/// Define the format constants
enum
{
  tile_format_jpeg = 0,
  tile_format_jpeg2000 = 1,
  tile_format_png = 2,
  tile_format_bmp = 3
};

/// VSF file content description
struct __attribute__((packed)) _index_file_content
{
  char header[30];                   // Should equal VMSCOPE_PRODUCT_HEADER
  uint8_t level_count;               // Number of pyramid layers
  uint8_t r;                         // RGB Background color R - channel
  uint8_t g;                         // RGB Background color G - channel
  uint8_t b;                         // RGB Background color B - channel
  int32_t size_x;                    // Total image width
  int32_t size_y;                    // Total image height
  int32_t resolution_x;              // X DPI Resolution
  int32_t resolution_y;              // Y DPI Resolution
  uint8_t format;                    // tile image format (see ETileFormat)
  uint8_t quality;                   // compression quality
  int32_t tile_size_x;               // single tile width
  int32_t tile_size_y;               // single tile height
  int32_t lowest_focal_plane_index;  // layer index of the lowest focal plane
  int32_t highest_focal_plane_index; // layer index of the highest focal plane
  float z_range;                     // distance between lowest and highest focal plane in micro meters
  uint8_t major_version;             // Product format major version
  uint8_t minor_version;             // Product format minor version
};
typedef struct _index_file_content index_file_content;

struct level_tile_data
{
  uint64_t offset;
  uint64_t size;
  uint64_t width;
  uint64_t height;
};

/**
 * Openslide management and data structs
 */
struct level
{
  struct _openslide_level base;
  struct _openslide_grid *grid;
  const char *filename;
  struct level_tile_data *tiles;

  int8_t layer;
  int64_t tiles_across;
  int64_t tiles_down;
};

struct vsf_ops_data
{
  index_file_content index_file_content;
  // struct vsf_slide_meta meta_data;
};

/*
 * Reusable helper functions
 */

/**
 * Get the number of digits in this number
 * @param number Number to determine the number of characters required to print this number
 * @return Number of characters required
 */
static inline uint8_t NUM_DIGITS(int32_t number)
{
  uint8_t result = (number < 0) ? 1 : 0;
  while (number != 0)
  {
    number /= 10;
    result++;
  }

  return result;
}

/**
 * Try to read the specified amount of bytes from the given file
 * @param fd File to read from
 * @param buffer Buffer to write data to - you have to ensure an appropriate buffer size
 * @param size Number of bytes to read
 * @return Number of bytes read
 */
static size_t _read_bytes(FILE *fd, void *buffer, size_t size)
{
  size_t bytes_read_total = 0;
  size_t blocks_read = 0;
  size_t bytes_remaining = size;
  size_t block_size = MIN(bytes_remaining, 4096);
  while ((blocks_read = fread(((char *)buffer) + bytes_read_total, block_size, bytes_remaining / block_size, fd)) > 0)
  {
    bytes_read_total += blocks_read * block_size;
    bytes_remaining -= blocks_read * block_size;
    block_size = MIN(bytes_remaining, 4096);

    if (bytes_remaining == 0)
    {
      break;
    }
  }

  return bytes_read_total;
}

/**
 * Try to read the specified amount of bytes from the given file
 * @param fd File to read from
 * @param size Number of bytes to read
 * @param offset Position in the file to seek to
 * @param err Error processor
 * @return Buffer with fetched data or NULL if an error occurred
 */
static char *_read_data(FILE *fd, size_t size, size_t offset, GError **err)
{
  // Allocate a buffer for input data
  char *buffer = g_slice_alloc(size);
  if (!buffer)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to create result buffer");
    return NULL;
  }

  // Fetch input data
  fseeko(fd, offset, SEEK_SET);
  if (_read_bytes(fd, buffer, size) != size)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to read data from file");
    return NULL;
  }

  return buffer;
}

/**
 * Return the content of the given node as NSString
 * @param node Node to extract content from
 * @return result string (you are responsile for freeing the string memory)
 */
static inline char *_create_string_from_node_content(const xmlNodePtr node)
{
  if (node == NULL)
  {
    return g_strdup("");
  }
  else
  {
    xmlChar *content = xmlNodeGetContent(node);
    if (content != NULL)
    {
      char *result = g_strdup((char *)content);
      xmlFree(content);

      return result;
    }
    else
    {
      return g_strdup("");
    }
  }
}

/**
 * Search the XML hierarchy for a single node with the given name
 * @param name Name to match
 * @param ns optional namespace the name must be in - to ignore namespace matching pass nil
 * @param root Node to use as root for searching
 * @return found node or NULL if nothing was found
 */
static inline xmlNodePtr _find_node_with_name(const char *name, const char *ns, const xmlNodePtr root)
{
  // You are only allowed to search element nodes
  if (root->type != XML_ELEMENT_NODE)
  {
    return NULL;
  }

  // Check the current node
  const xmlChar *element_name = root->name;
  const xmlChar *element_namespace = NULL;
  if (root->ns != NULL && root->ns->prefix != NULL)
  {
    element_namespace = root->ns->prefix;
  }

  if (xmlStrcmp(BAD_CAST name, element_name) == 0 && (ns == NULL || xmlStrcmp(BAD_CAST ns, element_namespace) == 0))
  {
    return root;
  }

  // Search children
  xmlNodePtr child = root->children;
  while (child != NULL)
  {
    xmlNodePtr result = _find_node_with_name(name, ns, child);
    if (result != NULL)
    {
      return result;
    }
    else
    {
      child = child->next;
    }
  }

  // Nothing found
  return NULL;
}

/**
 * Search the XML hierarchy for all nodes with the given name
 * @param name Name to match
 * @param ns optional namespace the name must be in - to ignore namespace matching pass nil
 * @param root Node to use as root for searching
 * @param result array where  found nodes should be written to
 */
static inline void _find_nodes_with_name(const char *name, const char *ns, const xmlNodePtr root, GPtrArray *result)
{
  // You are only allowed to search element nodes
  if (root->type != XML_ELEMENT_NODE)
  {
    return;
  }

  // Check the current node
  const xmlChar *element_name = root->name;
  const xmlChar *element_namespace = NULL;
  if (root->ns != NULL && root->ns->prefix != NULL)
  {
    element_namespace = root->ns->prefix;
  }

  if (xmlStrcmp(BAD_CAST name, element_name) == 0 && (ns == NULL || xmlStrcmp(BAD_CAST ns, element_namespace) == 0))
  {
    g_ptr_array_add(result, root);
  }

  // Search children
  xmlNodePtr child = root->children;
  while (child != NULL)
  {
    _find_nodes_with_name(name, ns, child, result);
    child = child->next;
  }
}

static inline void _color_from_string(const char *string, uint8_t *r, uint8_t *g, uint8_t *b)
{
  if (strcmp(string, "Black") == 0)
  {
    *r = 0;
    *g = 0;
    *b = 0;
  }
  else if (strcmp(string, "Blue") == 0)
  {
    *r = 0;
    *g = 0;
    *b = 255;
  }
  else if (strcmp(string, "Green") == 0)
  {
    *r = 0;
    *g = 255;
    *b = 0;
  }
  else if (strcmp(string, "Red") == 0)
  {
    *r = 255;
    *g = 0;
    *b = 0;
  }
  else if (strcmp(string, "Cyan") == 0)
  {
    *r = 0;
    *g = 255;
    *b = 255;
  }
  else if (strcmp(string, "Yellow") == 0)
  {
    *r = 255;
    *g = 255;
    *b = 0;
  }
  else if (strcmp(string, "Magenta") == 0)
  {
    *r = 255;
    *g = 0;
    *b = 255;
  }
  else if (strcmp(string, "White") == 0)
  {
    *r = 255;
    *g = 255;
    *b = 255;
  }
  else if (string[0] == '#')
  {
    uint32_t sr, sg, sb;
    sscanf(string, "#%02X%02X%02X", &sr, &sg, &sb);
    *r = sr;
    *g = sg;
    *b = sb;
  }
}

/*
 * Format specific data retrival operations
 */

/**
 * Read the index file from the given file pointer
 * This will consider the version of the file as specified in content
 * @param fd File to read from
 * @param content Index file information structure where data will be written to
 *                - an appropriate version must be set here when invoking this function
 * @param err Error container
 * @return true if reading the index file succeeded, else false
 */
static bool _read_index_file_content(FILE *fd, index_file_content *content, GError **err)
{
  /// Setup default values in the metadata struct
  content->r = 255;
  content->g = 255;
  content->b = 255;
  content->lowest_focal_plane_index = 0;
  content->highest_focal_plane_index = 0;
  content->quality = 0;
  content->format = tile_format_jpeg;
  content->z_range = 0;
  content->level_count = 9;
  content->resolution_x = 0;
  content->resolution_y = 0;

  if (content->major_version == 1)
  {
    off_t seek = 0;
    switch (content->minor_version)
    {
    case 0:
    {
      seek = 9;
      break;
    }
    case 1:
    {
      seek = 13;
      break;
    }
    case 2:
    {
      seek = 25;
      break;
    }
    default:
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unsupported product version");
      return false;
    }
    }

    // Seek to content start
    fseeko(fd, seek, SEEK_SET);
    if (_read_bytes(fd, &(content->size_x), sizeof(content->size_x)) != sizeof(content->size_x) ||
        _read_bytes(fd, &(content->size_y), sizeof(content->size_y)) != sizeof(content->size_y) ||
        _read_bytes(fd, &(content->tile_size_x), sizeof(content->tile_size_x)) != sizeof(content->tile_size_x) ||
        _read_bytes(fd, &(content->tile_size_y), sizeof(content->tile_size_y)) != sizeof(content->tile_size_y))
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed parsing header data");
      return false;
    }
  }
  else if (content->major_version == 2)
  {
    // Depending on the minor version is the available header information
    fseeko(fd, 0, SEEK_SET);
    size_t header_size = (content->minor_version == 0) ? 60 : 72;
    if (_read_bytes(fd, content, header_size) != header_size)
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed parsing header data");
      return false;
    }
  }

  return true;
}

/**
 * Read the VSF index file content
 * @param filename Filename to read the content from
 * @param content Result data structure where the fetched data will be pushed to
 * @param err Error pointer where any problems will be written to
 * @return true if the read was successfull else false
 */
static bool _read_index_file(const char *filename, index_file_content *content, GError **err)
{
  index_file_content result;

  // Validate file ending
  uint32_t length = strlen(filename);
  uint8_t extension_length = strlen(INDEX_FILE_EXTENSION);
  if (length <= extension_length)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Inappropriate filename");
    return false;
  }

  /// Validate extension
  char extension[extension_length];
  memcpy(extension, filename + (length - extension_length), extension_length * sizeof(char));
  uint8_t index = 0;
  for (index = 0; index < extension_length; index++)
  {
    // Lowercase the letter if required
    char s = extension[index] + ((extension[index] >= 'A' && extension[index] <= 'Z') ? 32 : 0);

    // Match with extension
    if (s != INDEX_FILE_EXTENSION[index])
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Inappropriate filename extension");
      return false;
    }
  }

  // Try to open the given filename
  FILE *vsf_file = _openslide_fopen(filename, "r", err);
  if (vsf_file == NULL)
  {
    return false;
  }

  // Try to read the header from the VSF file
  if (_read_bytes(vsf_file, result.header, 6) != 6)
  {
    fclose(vsf_file);
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to read product version");
    return false;
  }

  // Validate the product version
  boolean valid_product_version = false;
  if (result.header[1] == '1')
  {
    result.major_version = 1;

    switch (result.header[3])
    {
    case '0':
    case '1':
    case '2':
    {
      result.minor_version = (result.header[3] - '0');
      valid_product_version = true;
    }
    }
  }
  else if (result.header[3] >= '2' && result.header[5] >= '0' && result.header[5] <= '9')
  {
    result.major_version = 2;
    result.minor_version = (result.header[5] - '0');
    valid_product_version = true;
  }

  if (!valid_product_version)
  {
    fclose(vsf_file);
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to read product version");
    return false;
  }

  // Read the header
  boolean success = _read_index_file_content(vsf_file, &result, err);
  fclose(vsf_file);
  *content = result;
  return success;
}

/**
 * Build the filename with the additional extension
 * @param filename VSF filename
 * @param extension extension to append
 * @return Filename reflecting the given properties
 */
static char *_create_file_name_with_extension(const char *filename, const char *extension)
{
  uint32_t length = strlen(filename);
  uint8_t extension_length = strlen(INDEX_FILE_EXTENSION);

  // result length is the base filename length + the extension length
  uint32_t result_length = length - extension_length + strlen(extension);
  char *result = g_malloc0(result_length);
  memcpy(result, filename, (length - extension_length) * sizeof(char));
  memcpy(result + (length - extension_length), extension, strlen(extension) * sizeof(char));

  return result;
}

/**
 * Build the filename for the additional image source file based on the given
 * index file, image layer and focal plane index
 * @param file_info Index file information
 * @param filename VSF filename
 * @param layer Layer to get the filename for
 * @param focal_plane_index Focal plane index
 * @return Filename reflecting the given properties
 */
static char *_create_file_name_for_layer(index_file_content *file_info, const char *filename, uint8_t layer, int32_t focal_plane_index)
{
  // result length is at least the base filename length + the .img extension + -levelXX + terminating string
  uint32_t extension_length = strlen(IMAGE_FILE_EXTENSION) + strlen("-levelXX") + 2;
  if (focal_plane_index != 0) {
    // Focal plane id is infixed with + or - and the two digit id
    extension_length += 3;
  }

  char extension[extension_length];
  if (file_info->major_version == 1)
  {
    sprintf(extension, "-level%1i%s", layer, IMAGE_FILE_EXTENSION);
  }
  else if (focal_plane_index == 0)
  {
    sprintf(extension, "-level%02i%s", layer, IMAGE_FILE_EXTENSION);
  }
  else
  {
    sprintf(extension, "-level%02i%+02d%s", layer, focal_plane_index, IMAGE_FILE_EXTENSION);
  }

  return _create_file_name_with_extension(filename, extension);
}

/**
 * Check if the additional image source file is present and can be opened
 * @param file_info Index file information
 * @param filename filename VSF filename
 * @param layer Layer of the slide to open
 * @param focal_plane_index Focal plane index
 * @param err Error pointer where any problems will be written to
 * @return true if the file exists and is readable else false
 */
static bool _has_file_name_for_layer(index_file_content *file_info, const char *filename, uint8_t layer, int32_t focal_plane_index, GError **err)
{
  // Build the filename for the requested object
  char *image_filename = _create_file_name_for_layer(file_info, filename, layer, focal_plane_index);
  bool result = false;

  // Try to open the requested file
  FILE *file_pointer = _openslide_fopen(image_filename, "r", err);
  if (file_pointer != NULL)
  {
    result = true;
    fclose(file_pointer);
  }

  free(image_filename);
  return result;
}

/**
 * Fetch the tile information (location and size in the source file) from a version 1 file
 * @param minor_version Minor version information
 * @param fd File to read from
 * @param layer Layer the tile is located in
 * @param tile_index Index number of the tile in the layer
 * @param offset Result pointer where the tile offset is written to
 * @param size Result pointer where the tile size is written to
 * @param err Error pointer where any problems will be written to
 * @return true if the data could be extracted
 */
static bool _get_tile_infomation_version1(uint8_t minor_version, FILE *fd, uint8_t layer, uint32_t tile_index, uint64_t *offset, uint64_t *size, GError **err)
{
  off_t seek = 0;
  uint32_t tiles_x = 0, tiles_y = 0;
  uint8_t tile_record_size = 0, level_record_offset = 0, offset_size = 0;
  *offset = 0;
  *size = 0;

  switch (minor_version)
  {
  case 0:
  {
    seek = 25;
    tile_record_size = 12;
    level_record_offset = 16;
    offset_size = 4;
    break;
  }
  case 1:
  {
    seek = 29;
    tile_record_size = 16;
    level_record_offset = 16;
    offset_size = 8;
    break;
  }
  case 2:
  {
    seek = 41;
    tile_record_size = 16;
    level_record_offset = 28;
    offset_size = 8;
    break;
  }
  default:
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unsupported product version");
    return false;
  }
  }

  // Seek to the desired position and read the amount of tiles
  fseeko(fd, seek, SEEK_SET);
  if (_read_bytes(fd, &tiles_x, sizeof(tiles_x) != sizeof(tiles_x)) ||
      _read_bytes(fd, &tiles_y, sizeof(tiles_y) != sizeof(tiles_y)))
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed reading tile layout");
    return false;
  }

  // Proceed the file pointer to the desired level
  for (int i = 0; i < layer; i++)
  {
    fseeko(fd, tiles_x * tiles_y * tile_record_size + level_record_offset, SEEK_CUR);
  }

  // Validate data
  if (tiles_x * tiles_y <= tile_index)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Tile index is invalid - Number of tiles in file: %d - Index requested: %d", tiles_x * tiles_y, tile_index);
    return false;
  }

  // Skip to the position of the tile
  fseeko(fd, tile_index * tile_record_size, SEEK_CUR);
  if (_read_bytes(fd, offset, offset_size != sizeof(offset_size)))
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed reading tile offset");
    return false;
  }
  *offset = (*offset) >> ((sizeof(uint64_t) - offset_size) * 8);

  if (_read_bytes(fd, size, sizeof(uint32_t) != sizeof(uint32_t)))
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed reading tile size");
    return false;
  }
  *size = (*offset) >> 32;

  return true;
}

/**
 * Fetch the tile information (location and size in the source file) from a version 2 file
 * @param minor_version Minor version information
 * @param fd File to read from
 * @param tile_index Index number of the tile in the layer
 * @param offset Result pointer where the tile offset is written to
 * @param size Result pointer where the tile size is written to
 * @param err Error pointer where any problems will be written to
 * @return true if the data could be extracted
 */
static bool _get_tile_infomation_version2(uint8_t minor_version G_GNUC_UNUSED, FILE *fd, uint32_t tile_index, uint64_t *offset, uint64_t *size, GError **err)
{
  // Validate data
  fseeko(fd, sizeof(int64_t), SEEK_SET);
  uint64_t tile_count;
  if (_read_bytes(fd, &tile_count, sizeof(uint64_t)) != sizeof(uint64_t))
  {
    g_prefix_error(err, "Failed to read tile count from");
    return false;
  }

  if (tile_count <= tile_index)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Tile index is invalid - Number of tiles in file: %lu - Index requested: %d", tile_count, tile_index);
    return false;
  }

  // Read tile size
  fseeko(fd, tile_index * sizeof(uint64_t), SEEK_CUR);
  if (_read_bytes(fd, offset, sizeof(uint64_t)) != sizeof(uint64_t))
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to read tile offset");
    return false;
  }

  // Determine the tile size
  uint64_t next_tile_offset = 0;
  if (tile_index != tile_count - 1)
  {
    // Determine the location of the next tile in the file
    if (_read_bytes(fd, &next_tile_offset, sizeof(uint64_t)) != sizeof(uint64_t))
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to read follow up tile offset");
      return false;
    }
  }
  else
  {
    fseeko(fd, 0, SEEK_END);
    next_tile_offset = ftello(fd);
  }

  *size = next_tile_offset - *offset;
  return true;
}

/**
 * Fetch the tile information (location and size in the source file)
 * @param file_info Index file information
 * @param filename Filename where the tile data is stored in
 * @param layer Layer the tile is located in
 * @param tile_index Index number of the tile in the layer
 * @param offset Result pointer where the tile offset is written to
 * @param size Result pointer where the tile size is written to
 * @param err Error pointer where any problems will be written to
 * @return true if the data could be extracted
 */
static bool _get_tile_file_location(index_file_content *file_info, const char *filename, uint8_t layer, int32_t tile_index, uint64_t *offset, uint64_t *size, GError **err)
{
  // Try to open the desired file
  FILE *fd = _openslide_fopen(filename, "r", err);
  if (!fd)
  {
    g_prefix_error(err, "Can't open associated image %s: ", filename);
    return false;
  }

  bool success = false;
  if (file_info->major_version == 1)
  {
    success = _get_tile_infomation_version1(file_info->minor_version, fd, layer, tile_index, offset, size, err);
  }
  else if (file_info->major_version == 2)
  {
    success = _get_tile_infomation_version2(file_info->minor_version, fd, tile_index, offset, size, err);
  }
  else
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unsupported product version");
  }

  fclose(fd);
  return success;
}

/**
 * Get the dimensions of tile
 * @param osr Openstack whole slide pointer
 * @param level_data information about the requested resolution level
 * @param tile_index index of the tile in the given level
 * @param tile_col Column of the tile on the grid
 * @param tile_row Row of the tile on the grid
 * @param width Result pointer to write the tile width to
 * @param height Result pointer to write the tile height to
 * @param err Error pointer to write any errors to
 * @return true if the process was successfull else false
 */
static bool _get_tile_dimension(openslide_t *osr, struct level *level_data, int64_t tile_index, int64_t tile_col, int64_t tile_row, uint64_t *width, uint64_t *height, GError **err)
{
  // Determine the format of the tile
  struct vsf_ops_data *data = (struct vsf_ops_data *)osr->data;
  if (data->index_file_content.format == tile_format_jpeg)
  {
    int32_t tw, th;
    bool success = _openslide_jpeg_read_dimensions(level_data->filename, level_data->tiles[tile_index].offset, &tw, &th, err);
    *width = tw;
    *height = th;

    return success;
  }
  else
  {
    // Calculate the dimensions
    *width = MIN(level_data->base.tile_w, level_data->base.w - tile_col * level_data->base.tile_w);
    *height = MIN(level_data->base.tile_h, level_data->base.h - tile_row * level_data->base.tile_h);
    return true;
  }
}

/**
 * Decode a single tile in version 1 data format
 * @param layer_filename File to load data from
 * @param offset Offset in layer_filename where the tile data begins
 * @param size Amount of bytes in the input file from offset that belong to the tile
 * @param width Tile width
 * @param height Tile height
 * @param dest Destination buffer to write the result to
 * @param err Error handling pointer
 * @return true if the read was successfull, else no
 */
static bool _get_tile_data_version1(const char *layer_filename, uint64_t offset, uint64_t size, int64_t width, int64_t height, uint32_t *dest, GError **err)
{
  FILE *fd = _openslide_fopen(layer_filename, "r", err);
  if (!fd)
  {
    g_prefix_error(err, "Unable to open source file %s", layer_filename);
    return false;
  }

  /// Version 1 tiles are always jpeg encoded
  unsigned char jpegHeader[] = {0xff, 0xd8, 0xff, 0xe0, 0x00, 0x10, 0x4a, 0x46, 0x49, 0x46};
  unsigned char *buffer = g_slice_alloc(sizeof(jpegHeader) + size);

  fseeko(fd, offset, SEEK_SET);
  if (_read_bytes(fd, buffer + sizeof(jpegHeader), size) != size)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unable to read required amount of data");
    fclose(fd);
    return false;
  }

  // We finished our file transaction
  fclose(fd);

  bool decode_result = _openslide_jpeg_decode_buffer(buffer, size, dest, width, height, err);
  if (!decode_result)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to decode data");
  }

  g_slice_free1(sizeof(jpegHeader) + size, buffer);
  return decode_result;
}

/**
 * Decode a single tile in version 2 data format
 * @param index_file_content VSF meta data information
 * @param layer_filename File to load data from
 * @param offset Offset in layer_filename where the tile data begins
 * @param size Amount of bytes in the input file from offset that belong to the tile
 * @param width Tile width
 * @param height Tile height
 * @param dest Destination buffer to write the result to
 * @param err Error handling pointer
 * @return true if the read was successfull, else no
 */
static bool _get_tile_data_version2(index_file_content *index_file_content, const char *layer_filename, uint64_t offset, uint64_t size, int64_t width, int64_t height, uint32_t *dest, GError **err)
{
  FILE *fd = _openslide_fopen(layer_filename, "r", err);
  if (!fd)
  {
    g_prefix_error(err, "Unable to open source file %s", layer_filename);
    return false;
  }

  // Allocate a buffer for input data
  boolean needs_prefetch = (index_file_content->format == tile_format_jpeg || index_file_content->format == tile_format_jpeg2000);
  char *buffer = NULL;

  // JPEG and JPEG2000 can read data from buffers
  if (needs_prefetch == true)
  {
    buffer = _read_data(fd, size, offset, err);
    if (!buffer)
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unable to read required amount of data");
      fclose(fd);
      return false;
    }
  }

  // Decode input data based on the tile format
  boolean decode_result = false;
  switch (index_file_content->format)
  {
  case tile_format_jpeg:
  {
    decode_result = _openslide_jpeg_decode_buffer(buffer, size, dest, width, height, err);
    break;
  }
  case tile_format_jpeg2000:
  {
    decode_result = _openslide_jp2k_decode_buffer(dest, width, height, buffer, size, OPENSLIDE_JP2K_RGB, err);
    break;
  }
  case tile_format_png:
  {
    decode_result = _openslide_png_read(layer_filename, offset, dest, width, height, err);
    break;
  }
  case tile_format_bmp:
  {
    decode_result = _openslide_gdkpixbuf_read("bmp", layer_filename, offset, size, dest, width, height, err);
    break;
  }
  default:
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unknown tile data format");
    decode_result = false;
  }
  }

  // Free buffer
  if (buffer)
  {
    g_slice_free1(size, buffer);
  }

  fclose(fd);
  return decode_result;
}

static bool _get_tile_data(openslide_t *osr, struct level *level_data, int64_t tile_index, int64_t tile_width, int64_t tile_height, uint32_t *dest, GError **err)
{
  // Try to determine basic properties of the tile
  struct vsf_ops_data *data = (struct vsf_ops_data *)osr->data;
  const char *filename = g_hash_table_lookup(osr->properties, OPENSLIDE_PROPERTY_VSF_FILENAME);
  char *layer_filename = _create_file_name_for_layer(&data->index_file_content, filename, level_data->layer, 0);

  // Determine the offset and size of the tile
  uint64_t offset = 0, size = 0;
  if (!_get_tile_file_location(&data->index_file_content, layer_filename, level_data->layer, tile_index, &offset, &size, err))
  {
    g_free(layer_filename);
    return false;
  }

  bool result = false;
  if (size > 0)
  {
    if (data->index_file_content.major_version == 1)
    {
      result = _get_tile_data_version1(layer_filename, offset, size, tile_width, tile_height, dest, err);
    }
    else if (data->index_file_content.major_version == 2)
    {
      result = _get_tile_data_version2(&data->index_file_content, layer_filename, offset, size, tile_width, tile_height, dest, err);
    }
    else
    {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Unsupported product version");
    }
  }

  return result;
}

/**
 * Free all associated resources
 * @param osr
 */
static void destroy(openslide_t *osr)
{
  for (uint8_t layer = 0; layer < osr->level_count; layer++)
  {
    struct level *level_data = (struct level *)osr->levels[layer];
    g_free((void *)level_data->filename);
    g_slice_free1(sizeof(struct level_tile_data) * level_data->tiles_across * level_data->tiles_down, level_data->tiles);
  }

  struct vsf_ops_data *data = (struct vsf_ops_data *)osr->data;

  if (data)
  {
    g_slice_free(struct vsf_ops_data, data);
    osr->data = NULL;
  }
}

/**
 * Try to detect the VSF format in the given filename
 * @param filename VSF index file to open
 * @param tl TIFF-Like prefetched file content (must be NULL)
 * @param err Error pointer where any problems will be written to
 * @return true if this seams to be the VSF format else false
 */
static bool vsf_detect(const char *filename, struct _openslide_tifflike *tl, GError **err)
{
  // reject TIFFs
  if (tl)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Is a TIFF file");
    return false;
  }

  // Try to open fetch the metadata
  index_file_content content;
  if (!_read_index_file(filename, &content, err))
  {
    return false;
  }

  // Validate the presence of all required files
  for (uint8_t level = 0; level < content.level_count; level++)
  {
    for (int32_t focal_plane = content.lowest_focal_plane_index; focal_plane < content.highest_focal_plane_index; focal_plane++)
    {
      if (_has_file_name_for_layer(&content, filename, level, focal_plane, err) == false)
      {
        g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Missing image chunk file");
        return false;
      }
    }
  }

  return true;
}

/**
 * Paint the specified region to the given cairo surface
 * @param osr openslide information container
 * @param cr cairo surface to draw to
 * @param x x location to start drawing from
 * @param y y location to start drawing from
 * @param level Desired pyramid layer
 * @param w width to draw
 * @param h height to draw
 * @param err Error reporting container
 * @return true if drawing was successful else false
 */
static bool paint_region(openslide_t *osr G_GNUC_UNUSED, cairo_t *cr, int64_t x, int64_t y, struct _openslide_level *level, int32_t w, int32_t h, GError **err)
{
  struct level *l = (struct level *)level;
  bool success = _openslide_grid_paint_region(l->grid, cr, NULL, x / l->base.downsample, y / l->base.downsample, level, w, h, err);

  return success;
}

/**
 * Fetch the tile data
 * @param osr openslide information container
 * @param cr cairo surface to draw to
 * @param level Desired pyramid layer the tile is located in
 * @param tile_col X coordinate of the tile
 * @param tile_row Y coordinate of the tile
 * @param arg Additional arguments (unused)
 * @param err Error reporting container
 * @return true if data retrieval was successful else false
 */
static bool read_tile(openslide_t *osr, cairo_t *cr, struct _openslide_level *level, int64_t tile_col, int64_t tile_row, void *arg G_GNUC_UNUSED, GError **err)
{
  // Read the tile data
  struct level *level_data = (struct level *)level;
  int64_t tile_index = tile_row * level_data->tiles_across + tile_col;

  // Try to load data from cache
  struct _openslide_cache_entry *cache_entry;
  uint32_t *tiledata = _openslide_cache_get(osr->cache, level, tile_col, tile_row, &cache_entry);

  // Deterine tile dimensions (if not availble)
  struct level_tile_data *tile_data = level_data->tiles + tile_index;

  // No cached record was found, we need to load the data from file
  if (!tiledata)
  {
    struct vsf_ops_data *data = (struct vsf_ops_data *)osr->data;

    // Determine the offset and size of the tile
    if (!_get_tile_file_location(&data->index_file_content, level_data->filename, level_data->layer, tile_index, &tile_data->offset, &tile_data->size, err))
    {
      return false;
    }

    // Fetch tile dimensions
    if (!_get_tile_dimension(osr, level_data, tile_index, tile_col, tile_row, &tile_data->width, &tile_data->height, err))
    {
      return false;
    }

    // Fetch file properties
    int64_t buffer_size = tile_data->width * tile_data->height * 4;
    tiledata = g_slice_alloc(buffer_size);
    _get_tile_data(osr, level_data, tile_index, tile_data->width, tile_data->height, tiledata, err);
    if (openslide_get_error(osr) != NULL)
    {
      g_slice_free1(buffer_size, tiledata);
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Failed to read tile data");
      return false;
    }

    // put it in the cache
    _openslide_cache_put(osr->cache, level, tile_col, tile_row, tiledata, buffer_size, &cache_entry);
  }

  // draw it
  cairo_surface_t *surface = cairo_image_surface_create_for_data((unsigned char *)tiledata, CAIRO_FORMAT_ARGB32, tile_data->width, tile_data->height, tile_data->width * 4);
  cairo_set_source_surface(cr, surface, 0, 0);
  cairo_surface_destroy(surface);
  cairo_paint(cr);

  // done with the cache entry, release it
  _openslide_cache_entry_unref(cache_entry);

  return true;
}

/**
 * Sorting callback for pyramid layers
 * @param a Layer to compare
 * @param b Layer to compare
 * @return -1 if a is wider than b, 0 if they are of equal width or 1 if b is wider than a
 */
static int width_compare(gconstpointer a, gconstpointer b)
{
  const struct level *la = *(const struct level **)a;
  const struct level *lb = *(const struct level **)b;

  if (la->base.w > lb->base.w)
  {
    return -1;
  }
  else if (la->base.w == lb->base.w)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

/* the function pointer structure for backends */
static const struct _openslide_ops vsf_ops = {
    .paint_region = paint_region,
    .destroy = destroy,
};

static bool vsf_open(openslide_t *osr, const char *filename, struct _openslide_tifflike *tl, struct _openslide_hash *quickhash1 G_GNUC_UNUSED, GError **err)
{
  // Read the meta file to obtain information
  // We can savely assume that everything is formatted correctly, as the vsf_detect already validated the index file
  struct vsf_ops_data *data = g_slice_new0(struct vsf_ops_data);

  if (!_read_index_file(filename, &data->index_file_content, err))
  {
    return false;
  }

  // Reject TIFFs
  if (tl)
  {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED, "Is a TIFF file");
    return false;
  }

  osr->level_count = data->index_file_content.level_count;

  // Build the levels array
  GPtrArray *level_array = g_ptr_array_new();
  for (uint8_t layer = 0; layer < data->index_file_content.level_count; layer++)
  {
    struct level *level_data = g_slice_new0(struct level);
    level_data->base.w = data->index_file_content.size_x >> layer;
    level_data->base.h = data->index_file_content.size_y >> layer;
    level_data->base.tile_w = data->index_file_content.tile_size_x;
    level_data->base.tile_h = data->index_file_content.tile_size_y;
    level_data->layer = layer;
    level_data->filename = _create_file_name_for_layer(&data->index_file_content, filename, layer, 0);
    level_data->tiles_across = (int)ceilf((float)level_data->base.w / (float)level_data->base.tile_w);
    level_data->tiles_down = (int)ceilf((float)level_data->base.h / (float)level_data->base.tile_h);
    level_data->tiles = g_slice_alloc0(sizeof(struct level_tile_data) * level_data->tiles_across * level_data->tiles_down);

    // create grid for our level
    level_data->grid = _openslide_grid_create_simple(osr, level_data->tiles_across, level_data->tiles_down, data->index_file_content.tile_size_x, data->index_file_content.tile_size_y, read_tile);
    g_ptr_array_add(level_array, level_data);
  }

  // sort tiled levels
  g_ptr_array_sort(level_array, width_compare);
  struct level *top_level = g_ptr_array_index(level_array, 0);

  osr->levels = (struct _openslide_level **)g_ptr_array_free(level_array, false);
  level_array = NULL;

  osr->data = data;

  // Setup properties
  g_hash_table_insert(osr->properties, g_strdup(OPENSLIDE_PROPERTY_NAME_COMMENT), g_strdup(data->index_file_content.header));
  g_hash_table_insert(osr->properties, g_strdup(OPENSLIDE_PROPERTY_VSF_FILENAME), g_strdup(filename));
  g_hash_table_insert(osr->properties, g_strdup(OPENSLIDE_PROPERTY_NAME_MPP_X), _openslide_format_double(1.0 / ((double)data->index_file_content.resolution_x / 25400.0)));
  g_hash_table_insert(osr->properties, g_strdup(OPENSLIDE_PROPERTY_NAME_MPP_Y), _openslide_format_double(1.0 / ((double)data->index_file_content.resolution_y / 25400.0)));

  _openslide_set_background_color_prop(osr, data->index_file_content.r, data->index_file_content.g, data->index_file_content.b);
  _openslide_set_bounds_props_from_grid(osr, top_level->grid);

  // Assign paint and destroy operations
  osr->ops = &vsf_ops;

  return true;
}

/// public exported methods

const struct _openslide_format _openslide_format_vsf = {
    .name = "vsf",
    .vendor = "vsf",
    .detect = vsf_detect,
    .open = vsf_open,
};
