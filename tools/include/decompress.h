/******************************************************************************\
* SOUBOR: decompress.h
* AUTOR: Drahomir Dlabaja (xdlaba02)
\******************************************************************************/

#ifndef DECOMPRESS_H
#define DECOMPRESS_H

#include <cstdint>
#include <vector>
#include <string>

using namespace std;

void print_usage(const char *argv0);

bool parse_args(int argc, char *argv[], const char *&input_file_name, const char *&output_file_mask);

bool savePPMs(const vector<uint8_t> &rgb_data, uint64_t width, uint64_t height, uint32_t color_image_count, uint64_t image_count, const string &output_file_mask);

#endif
