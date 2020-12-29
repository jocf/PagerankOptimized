#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>

#include "./lib/pagerank.h"

// This method continually returns 1 until the convergence threshold is reach when it will return 0
// This is the parallelized method that is more performant and is the method actually used when the program is run
double iterate_pagerank_optimized(list* plist, double **pageranks, double dampener, int npages, int swap_value, int n_threads){
	// if swap value is == 0, then pageranks[1] is the list to update and compare with pageranks[0]
	// if swap value is == 1, then pageranks[0] is the list to update and compare with pageranks[1]
	
	omp_set_num_threads(n_threads);
	
	register int write_index = 0;
	register int read_index = 1;
	
	if(swap_value == 0){
		write_index = 1;
		read_index = 0;
	}
	
	register struct node** node_list = malloc(sizeof(struct node*)*npages);
	register struct node* tmp_current_node = plist->head;
	
	// We will transfer to a dynamically allocated array of pointers, as using a linked list
	// with a multithreaded implementation is extremely problematic. 
	
	for(int i = 0; i < npages-1; i+=2){
		node_list[i] = tmp_current_node;
		tmp_current_node = tmp_current_node->next;
		node_list[i+1] = tmp_current_node;
		tmp_current_node = tmp_current_node->next;
	}
	
	if(npages%2 != 0){
		node_list[npages-1] = tmp_current_node;
	}
	
	
	#pragma omp parallel for schedule(auto)
	// Convert this to an array and access the correct member with the i value for indexing inside the loop
	for (int i = 0; i < npages; i++){
		register struct node* current_node = node_list[i];
		// Now that we have fetched the correct node, we can perform the calculation.
		
		register double new_value = 0.0;
		
		new_value += (1.0 - dampener)/npages;
		
		if(current_node->page->inlinks != NULL){
			
			register struct node* inlinks_node = current_node->page->inlinks->head;
			
			for(int j = 0; j < current_node->page->inlinks->length - 1; j+=2){
				//struct node* inlinks_node = current_node->page->inlinks->head;
				// fetch the inlink page 
				
				new_value += dampener * (pageranks[read_index][inlinks_node->page->index] / inlinks_node->page->noutlinks);
				inlinks_node = inlinks_node->next;
				new_value += dampener * (pageranks[read_index][inlinks_node->page->index] / inlinks_node->page->noutlinks);
				inlinks_node = inlinks_node->next;
				
				
			}
			if(current_node->page->inlinks->length %2 != 0){
				new_value += dampener * (pageranks[read_index][inlinks_node->page->index] / inlinks_node->page->noutlinks);
			}
			
		}
		
		pageranks[write_index][i] = new_value;
	
	}
	free(node_list);
	
	
	double convergence_value = 0.0;
	#pragma omp parallel
	{
		//private variables here
		register double priv_converge = 0.0;
		register double tmp_1 = 0.0;
		register double tmp_2 = 0.0;
		#pragma omp for schedule(auto)
		for(int i = 0; i < npages - 1; i +=2 ){
			tmp_1 += pageranks[write_index][i] - pageranks[read_index][i];
			tmp_1 = tmp_1 * tmp_1;
			tmp_2 += pageranks[write_index][i+1] - pageranks[read_index][i+1];
			tmp_2 = tmp_2 * tmp_2;
			
			priv_converge += tmp_2 + tmp_1;
			tmp_1 = 0.0;
			tmp_2 = 0.0;	
		}
		
		// Critical was faster and more consistent than atomic after thorough testing
		#pragma omp critical
		convergence_value += priv_converge;
	}
	if(npages % 2 != 0){
		// do the last element
		double tmp_1 = 0.0;
		tmp_1 += pageranks[write_index][npages-1] - pageranks[read_index][npages-1];
		tmp_1 = tmp_1 * tmp_1;
		convergence_value += tmp_1;
	}
	
	convergence_value = sqrt(convergence_value);
	
	if(convergence_value <= EPSILON){
		register struct node* tmp_node = plist->head;
		for(int i = 0; i < npages - 1; i+=2){
			printf("%s %.4lf\n", tmp_node->page->name,pageranks[write_index][i]);
			tmp_node = tmp_node->next;
			printf("%s %.4lf\n", tmp_node->page->name,pageranks[write_index][i+1]);
			tmp_node = tmp_node->next;
		}
		
		if(npages % 2 != 0){
			printf("%s %.4lf\n", tmp_node->page->name,pageranks[write_index][npages-1]);
		}
		
		return 0;
	}
	
	return 1;
}


// This method continually returns 1 until the convergence threshold is reach when it will return 0
// This is the serial, unoptimized method created for testing and comparative purposes

double iterate_pagerank(list* plist, double **pageranks, double dampener, int npages, int swap_value, int num_threads){
	// if swap value is == 0, then pageranks[1] is the list to update and compare with pageranks[0]
	// if swap value is == 1, then pageranks[0] is the list to update and compare with pageranks[1]
	
	register int write_index = 0;
	register int read_index = 1;
	
	if(swap_value == 0){
		write_index = 1;
		read_index = 0;
	}
	struct node* current_node = plist->head;
	for (int i = 0; i < npages; i++){
		
		// Now that we have fetched the correct node, we can perform the calculation.
		
		register double new_value = 0.0;
		
		new_value += (1.0 - dampener)/npages;
		
		if(current_node->page->inlinks != NULL){
		
			
			struct node* inlinks_node = current_node->page->inlinks->head;
			for(int j = 0; j < current_node->page->inlinks->length; j++){
				//struct node* inlinks_node = current_node->page->inlinks->head;
				// fetch the inlink page 
				

				new_value += dampener * (pageranks[read_index][inlinks_node->page->index] / inlinks_node->page->noutlinks);
				inlinks_node = inlinks_node->next;
				
			}
		}
		pageranks[write_index][i] = new_value;
		current_node = current_node->next;
		
	
	}
	
	register double convergence_value = 0.0;
	if(swap_value == 0){
		// pageranks[1] is updated value, pageranks[0] is old value
		// therefore, we do pageranks[1] - pageranks[0]
		for(int i = 0; i < npages; i++){
			register double tmp = 0.0;
			tmp += pageranks[1][i] - pageranks[0][i];
			tmp = tmp * tmp;
			convergence_value += tmp;
		}
	}
	else{
		// pageranks[0] is updated value, pageranks[1] is old value
		// therefore, we do pageranks[0] - pageranks[1]
		for(int i = 0; i < npages; i++){
			register double tmp = 0.0;
			tmp += pageranks[0][i] - pageranks[1][i];
			tmp = tmp * tmp;
			convergence_value += tmp;
		}
		
	}
	
	
	convergence_value = sqrt(convergence_value);
	//printf("%f\n",convergence_value);
	
	if(convergence_value <= EPSILON){
		//printf("Done! \n");
		struct node* tmp_node = plist->head;
		if(swap_value == 0){
			for(int i = 0; i < npages; i++){
				printf("%s %.4lf\n", tmp_node->page->name,pageranks[1][i]);
				tmp_node = tmp_node->next;
			}
			
		}else{
			for(int i = 0; i < npages; i++){
				printf("%s %.4lf\n", tmp_node->page->name,pageranks[0][i]);
				tmp_node = tmp_node->next;
			}
			
		}
		
		return 0;
	}
	
	return 1;
}

void pagerank(list* plist, int ncores, int npages, int nedges, double dampener) {
	// The first step will be to create the 2d array to store the pagerank data in
	// We will use a 2 dimensional array to store the previous, and the next pagerank data
	// I.e. to store t, and t+1, so we can compare between iterations
	
	// swap values to registers that will be read to continually
	
	register int npages_reg = npages;
	
	register double dampener_reg = dampener;
	
	// The number of threads, after a lot of testing, has been decided to be equal to the ncores
	// specified by the test file
	
	
	int n_threads = ncores * 1;
	/* TO BE COMMENTED OUT BEFORE SUBMISSION
	if(n_threads > 4){
		n_threads = 4;
	}
	-------------------------------------*/
	
	register double **pageranks = malloc(sizeof(double*)*2);
	pageranks[0] = malloc(sizeof(double)*npages);
	pageranks[1] = malloc(sizeof(double)*npages);
	
	// Next, we need to initialize all values to 1.0/npages
	omp_set_num_threads(n_threads);
	#pragma omp parallel for
	for(int i = 0; i < npages - 1; i+= 2){
		pageranks[0][i] = 1.0/npages;
		pageranks[1][i] = 1.0/npages;
		pageranks[0][i+1] = 1.0/npages;
		pageranks[1][i+1] = 1.0/npages;
	}
	
	if(npages % 2 != 0){
		pageranks[0][npages-1] = 1.0/npages;
		pageranks[1][npages-1] = 1.0/npages;
	}
	// Used for determining which index is the up to date (t+1) index in the pageranks array
	int swap_value = 0;
	
	
	// Iterate the pageranker method
	
	// If you would like to test the unoptimized method, please change iterate_pagerank_optimized to iterate_pagerank below \/\/\/
	while(iterate_pagerank_optimized(plist, pageranks, dampener_reg, npages_reg, swap_value, n_threads)){
		
		if(swap_value == 0){
			
			swap_value = 1;
		}
		else if(swap_value == 1){
			
			swap_value = 0;
		}
	}
	
	free(pageranks[0]);
	free(pageranks[1]);
	free(pageranks);
    
}

/*
######################################
### DO NOT MODIFY BELOW THIS POINT ###
######################################
*/

int main(void) {

    /*
    ######################################################
    ### DO NOT MODIFY THE MAIN FUNCTION OR HEADER FILE ###
    ######################################################
    */

    list* plist = NULL;

    double dampener;
    int ncores, npages, nedges;

    /* read the input then populate settings and the list of pages */
    read_input(&plist, &ncores, &npages, &nedges, &dampener);

    /* run pagerank and output the results */
    pagerank(plist, ncores, npages, nedges, dampener);

    /* clean up the memory used by the list of pages */
    page_list_destroy(plist);

    return 0;
}
