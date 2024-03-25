#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated

void heapify(int pData[], int N, int i)
{
 
    // Initialize largest as root
    int largest = i;
 
    // left = 2*i + 1
    int l = 2 * i + 1;
 
    // right = 2*i + 2
    int r = 2 * i + 2;
 
    // If left child is larger than root
    if (l < N && pData[l] > pData[largest])
        largest = l;
 
    // If right child is larger than largest
    // so far
    if (r < N && pData[r] > pData[largest])
        largest = r;
 
    // If largest is not root
    if (largest != i) {
        int temp = pData[i];
        pData[i] = pData[largest];
        pData[largest] = temp;

        heapify(pData, N, largest);
    }
}

void heapSort(int arr[], int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for (int i = n - 1; i > 0; i--)
    {
        int temp = arr[0];
        arr[0] = arr[i];
        arr [i] = temp;

        heapify(arr, i, 0);
    }
}

// implement merge sort
void merge(int pData[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    int L[n1], R[n2];

    for (i = 0; i < n1; i++)
        L[i] = pData[l + i];
    for (j = 0; j < n2; j++)
        R[j] = pData[m + 1 + j];

    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            pData[k] = L[i];
            i++;
        }
        else
        {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        pData[k] = L[i];
        i++;
        k++;
    }

    while (j < n2)
    {
        pData[k] = R[j];
        j++;
        k++;
    }
}
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
	if (l < r)
    {
        int m = l + (r - l) / 2;

        mergeSort(pData, l, m);
        mergeSort(pData, m + 1, r);

        merge(pData, l, m, r);
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{
	for (int i = 1; i < n; ++i)
    {
        int key = pData[i];
        int j = i - 1;
        while (j >= 0 && pData[j] > key)
        {
            pData[j + 1] = pData[j];
            j = j - 1;
        }
        pData[j + 1] = key;
    }
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
	for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < n - i - i; ++j)
        {
            if (pData[j] > pData[j + i])
            {
                int temp = pData[j];
                pData[j] = pData[j + i];
                pData[j + i] = temp;
            }
        }
    }
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
	for (int i = 0; i < n -1; ++i)
    {
        int minIndex = 1;
        for (int j = i +1; j < n; ++j)
        {
            if (pData[j] < pData[minIndex])
            {
                minIndex = j;
            }
        }

        int temp = pData[i];
        pData[i] = pData[minIndex];
        pData[minIndex] = temp;
    }
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	FILE* inFile = fopen(inputFileName,"r");
	int dataSz = 0;
	int i, n, *data;
	*ppData = NULL;
	
	if (inFile)
	{
		fscanf(inFile,"%d\n",&dataSz);
		*ppData = (int *)malloc(sizeof(int) * dataSz);
		// Implement parse data block
		if (*ppData == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(-1);
		}
		for (i=0;i<dataSz;++i)
		{
			fscanf(inFile, "%d ",&n);
			data = *ppData + i;
			*data = n;
		}

		fclose(inFile);
	}
	
	return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz );
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
