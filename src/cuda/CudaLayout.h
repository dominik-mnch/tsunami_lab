#ifndef CUDA_LAYOUT_H
#define CUDA_LAYOUT_H

// ROW-MAJOR:
// Index = row * stride + col
// stride = nCellsX + 2 (two ghost cells)
#define ROW_MAJOR(row, col, stride) ( (row) * (stride) + (col) )

#endif