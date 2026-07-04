// CudaLayout.h
#ifndef CUDA_LAYOUT_H
#define CUDA_LAYOUT_H

#define ROW_MAJOR(row, col, stride) ((row) * (stride) + (col))
#define COLUMN_MAJOR(row, col, height) ((col) * (height) + (row))

#endif 

