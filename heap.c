/**
 * 最小堆函数
 * 修改自开源项目：
 * https://github.com/armon/c-minheap-array
 * 作者：cangye@hotmail.com
 * 仅用于快速行进算法寻找最小值
 */

#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "../include/heap.h"
#include <stdio.h>
// 堆构建宏
#define LEFT_CHILD(i)   ((i<<1)+1)
#define RIGHT_CHILD(i)  ((i<<1)+2)
#define PARENT_ENTRY(i) ((i-1)>>1)
#define SWAP_ENTRIES(parent,child)  { \
                                      float temp1 = parent->key; \
                                      parent->key = child->key;          \
                                      child->key = temp1;                 \
                                      mincoord temp2 = parent->value;              \
                                      parent->value = child->value;      \
                                      child->value = temp2;               \
                                    }
#define GET_ENTRY(index,table) ((heap_entry*)(table+index))
// 比较函数
int compare_int_keys(register float key1_v, register float key2_v) {
    if (key1_v < key2_v)
        return -1;
    else if (key1_v == key2_v)
        return 0;
    else
        return 1;
}
void heap_create(heap* h, int initial_size, int (*comp_func)(float, float)){
    // 设置最小值大小
    int MIN_ENTRY_SIZE = 1024;
    // 检查初始大小
    if (initial_size <= MIN_ENTRY_SIZE)
        initial_size = MIN_ENTRY_SIZE;
    int MIN_MEM_SIZE = initial_size * sizeof(heap_entry);
    // 比较函数
    if (comp_func == NULL)
        comp_func = compare_int_keys;
    h->active_entries = 0;
    h->compare_func = comp_func;
    // 堆大小设置为0
    
    // 分配内存大小
    
    h->allocated_size = MIN_MEM_SIZE;
    h->minimum_size = h->allocated_size;
    // 分配内存
    //printf("分配内存c%d\n", h->allocated_size); 
    h->table = (void*)calloc(h->allocated_size, sizeof(heap_entry));
}

// 删除堆
void heap_destroy(heap* h) {
    // 错误检查
    assert(h != NULL);
    // 删除表
    free(h->table);
    // 删除其他内容
    h->active_entries=0;
    h->allocated_size = 0;
    h->table = NULL;
}

// 获取堆大小
int heap_size(heap* h) {
    return h->active_entries;
}

// 获取最小元素
int heap_min(heap* h, float* key, mincoord* value) {
    // 异常检测
    if (h->active_entries == 0)
        return 0;

    heap_entry* root = GET_ENTRY(0, h->table);

    // 获取最小值地址
    *key = root->key;
    *value = root->value;

    // 成功返回1 
    return 1;
}

// 插入新元素
void heap_push(heap *h, float key, mincoord value) {
    assert(h->table != NULL);
    // 现有内存大小
    int max_entries = h->allocated_size;
    // 防止频繁分配内存
    if (h->active_entries + 1 > max_entries) {
        // 如果新加元素大于所需内存
        int new_size = h->allocated_size * 2;
        // 分配新内存
        //printf("分配内存b%d\n", new_size); 
        heap_entry* new_table = (void*)calloc(new_size, sizeof(heap_entry));

        // 复制内存
        memcpy(new_table, h->table, h->allocated_size);
        
        // 释放内存
        free(h->table);
        // 更换内存
        h->table = new_table;
        h->allocated_size = new_size;
    }
    
    // 比较函数
    int (*cmp_func)(float, float) = h->compare_func;

    // 内存指针
    heap_entry* table = h->table;

    // 获取当前大小
    int current_index = h->active_entries;
    // 新元素到堆尾
    heap_entry* current = GET_ENTRY(current_index, table);

    // 循环
    int parent_index;
    heap_entry *parent;

    // 构建最小堆
    while (current_index > 0) {
        // 获取父节点索引
        parent_index = PARENT_ENTRY(current_index);

        // 获取父节点
        parent = GET_ENTRY(parent_index, table);
       
        // 比较值
        if (cmp_func(key, parent->key) < 0) {
            // 交换
            current->key = parent->key;
            current->value = parent->value;

            // 移动到父节点
            current_index = parent_index;
            current = parent;

        // 交换完成
        }   else
            break;
    }

    // 在当前节点插入
    current->key = key;
    current->value = value; 

    // 计算内存值
    h->active_entries++;
}


// 删除堆元素最小值
int heap_pop(heap* h, float* key, mincoord* value) {
    // 错误检查
    if (h->active_entries == 0)
        return 0;

    // 加载
    heap_entry* table = h->table;

    // 获取根节点
    int current_index = 0;
    heap_entry* current = GET_ENTRY(current_index, table);

    // 输出
    *key = current->key;
    *value = current->value;

    // 减少索引
    h->active_entries--;

    // 获取最后节点
    int entries = h->active_entries;
   
    // 判断是否是长度为1的堆，是的话不用处理
    if (h->active_entries > 0) {
        // 把最后节点赋予初始节点
        heap_entry* last = GET_ENTRY(entries,table);
        current->key = last->key;
        current->value = last->value;

        // 左子树和右子树
        heap_entry* left_child;
        heap_entry* right_child;

        // 比较函数
        int (*cmp_func)(float, float) = h->compare_func;

        // 左子树索引
        int left_child_index;
        // 迭代调整为小根堆
        while (left_child_index = LEFT_CHILD(current_index), left_child_index < entries) {
            // 左子树
            left_child = GET_ENTRY(left_child_index, table);
            if (left_child_index+1 < entries) {
                // 右子树
                right_child = GET_ENTRY((left_child_index+1), table);

                // 获取较小值
                if (cmp_func(left_child->key, right_child->key) <= 0) {

                    // 交换
                    if (cmp_func(current->key, left_child->key) == 1) {
                        SWAP_ENTRIES(current,left_child);
                        current_index = left_child_index;
                        current = left_child;

                    // 
                    } else
                        break;

                // 
                } else {

                    // 交换
                    if (cmp_func(current->key, right_child->key) == 1) {
                        SWAP_ENTRIES(current,right_child);
                        current_index = left_child_index+1;
                        current = right_child;

                    // 无
                    } else
                        break;

                }


            // 如果仅有左子树
            } else if (cmp_func(current->key, left_child->key) == 1) {
                SWAP_ENTRIES(current,left_child);
                current_index = left_child_index;
                current = left_child;

            // 无
            }  else
                break;

        }
    } 

    // 检查是否需要释放内存
    int used_size = entries * sizeof(heap_entry);

    // Allow one empty page, but not two
    if (h->allocated_size / 2 > used_size + 1 && h->allocated_size / 2 >= h->minimum_size) {
        // 计算需要页大小
        int new_size = h->allocated_size / 2;

        // 分配新内存
        //printf("分配内存a%d\n", new_size); 
        heap_entry* new_table = (void*)calloc(new_size, sizeof(heap_entry));

        // 拷贝
        memcpy(new_table, h->table, used_size);
        
        // 清除内存
        free(h->table);
        // 无
        h->table = new_table;
        h->allocated_size = new_size;
    }

    // 成功
    return 1;
}
