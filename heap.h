#ifndef MYHEAP_H
#define MYHEAP_H

// 构建坐标元素
typedef struct mincoord {
    int n1, n2, n3; 
} mincoord; 
typedef struct heap_entry {
    float key;   // 对比值
    mincoord value; // 取值
} heap_entry;

// 堆结构
typedef struct heap {
    int (*compare_func)(float, float); // 对比函数
    int minimum_size;   // 最小内存大小
    int active_entries;    // 求活跃块大小（堆大小）
    int allocated_size; // 内存大小
    heap_entry* table; // 堆内存结构
} heap;


/**
 * Creates a new heap
 * @param h Pointer to a heap structure that is initialized
 * @param initial_size What should the initial size of the heap be. If <= 0, then it will be set to the minimum
 * permissable size, of 1 page (512 entries on 32bit system with 4K pages).
 * @param comp_func A pointer to a function that can be used to compare the keys. If NULL, it will be set
 * to a function which treats keys as signed ints. This function must take two keys, given as pointers and return an int.
 * It should return -1 if key 1 is smaller, 0 if they are equal, and 1 if key 2 is smaller.
 */
void heap_create(heap* h, int initial_size, int (*comp_func)(float, float));

/**
 * Returns the size of the heap
 * @param h Pointer to a heap structure
 * @return The number of entries in the heap.
 */
int heap_size(heap* h);

/**
 * Inserts a new element into a heap.
 * @param h The heap to insert into
 * @param key The key of the new entry
 * @param value The value of the new entry
 */
void heap_push(heap* h, float key, mincoord value);

/**
 * Returns the element with the smallest key in the heap.
 * @param h Pointer to the heap structure
 * @param key A pointer to a pointer, to set to the minimum key
 * @param value Set to the value corresponding with the key
 * @return 1 if the minimum element exists and is set, 0 if there are no elements.
 */
int heap_min(heap* h, float* key, mincoord* value);

/**
 * Deletes the element with the smallest key from the heap.
 * @param h Pointer to the heap structure
 * @param key A pointer to a pointer, to set to the minimum key
 * @param valu Set to the value corresponding with the key
 * @return 1if the minimum element exists and is deleted, 0 if there are no elements.
 */
int heap_pop(heap* h, float* key, mincoord* value);


/**
 * Destroys and cleans up a heap.
 * @param h The heap to destroy.
 */
void heap_destroy(heap* h);
#endif