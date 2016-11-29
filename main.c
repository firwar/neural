
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NUM_INPUTS 2
#define NUM_LAYERS 2
#define NUM_NODES 6
#define SPEED 0.25;

//#define PRINTALL

int inputseq;
double averageerror;
double inputseq0;
double inputseq1;
double inputseq2;
double inputseq3;

// Initializes random function generator
int init_system() {
	time_t t;
	/* Intializes random number generator */
	srand((unsigned) time(&t));
}

// Get the input training set for this matrix
int getInputs(double* inputs, int num) {

	inputseq = rand()%4;
	if (inputseq == 0) {
		inputs[0] = 0;
    	inputs[1] = 0;
	} else if (inputseq == 1) {
		inputs[0] = 0;
    	inputs[1] = 1;
	} else if (inputseq == 2) {
		inputs[0] = 1;
    	inputs[1] = 0;
	} else {
		inputs[0] = 1;
    	inputs[1] = 1;
	}

	//inputs[0] = 1;
	//inputs[1] = 1;

	for (int i=0; i<num; i++)
	{

	}

	return 0;
}

// Take a matrix of row rows and cols columns and assign
// random values between 0 and 1 to them
int randomizeWeights(double* weight_matrix, int rows, int cols) {
	double generated_num;
   	int n = 1000;

	// Assign random values between 0 to 1 with 1/n steps
	for (int i=0; i<cols; i++) {
		for (int j=0; j<rows; j++) {
				generated_num = rand() % n;
				generated_num = generated_num/n;
				weight_matrix[j+i*rows] = generated_num;
				printf("%f\n", weight_matrix[j+i*rows]);
		}
	}

	return 0;
}

// Apply activation sigmoid function to each node
double applyActivation(double x) {

     double exp_value;
     double sigmoid_output;
     /*** Exponential calculation ***/
     exp_value = exp((double) -x);

     /*** Final sigmoid value ***/
     sigmoid_output = 1.0 / (1.0 + exp_value);

     return sigmoid_output;
}

// Apply activation sigmoid function to each node
double applyActivationDerivative(double x) {

     /*** Final sigmoid value ***/
     double sigmoid_deriv_output = x * (1.0 - x);

     return sigmoid_deriv_output;
}

// Compute the inner layer values including activation function
int computeInnerLayers(double *inputs, double *input_weight_matrix, double *inner_layer_matrix) {

	for (int i=0; i<NUM_NODES; i++) {
		for (int j=0; j<NUM_INPUTS; j++) {
			
			// Initialize inner layer matrix to 0
			if (j==0) {
				inner_layer_matrix[i] = 0;
			}
			
			//Cumulative dot product 
			inner_layer_matrix[i] += (inputs[j] * input_weight_matrix[i+j*NUM_NODES]);
			#ifdef PRINTALL
			printf("Input Weight matrix =  %f \n", input_weight_matrix[i+j*NUM_NODES]);
			#endif
		}		
	}

    for (int i=0; i<NUM_NODES; i++) {
		inner_layer_matrix[i] = applyActivation(inner_layer_matrix[i]);
	}

	for (int i=0; i<NUM_NODES; i++) {
		#ifdef PRINTALL
		printf("Inner layer matrix =  %f \n", inner_layer_matrix[i]);
		#endif
	}
	return 0;

}

// Compute outputs
double computeOutput(double *output_weight_matrix, double *inner_layer_matrix) {
	double output = 0;
	for (int i=0; i<NUM_NODES; i++) {
			output += output_weight_matrix[i] * inner_layer_matrix[i];
	}
	output = applyActivation(output);
	#ifdef PRINTALL
	printf("Output = %f\n", output);
	#endif
	return output;
}


// Return the output of an XOR which is the target
double getTarget(double *inputs) {
	if (inputs[0] == inputs[1]) {
		return 0.01;
	} else {
		return 1;
	}
}

// Calculate the output distance to travel along the sigmoid
double getDeltaOutputSum(double target, double output) {
	
	double delta_output_sum;
	
	double error = sqrtf((target - output)*(target-output));
	double error_prime = -(target - output);
	delta_output_sum = error_prime * applyActivationDerivative(output);

	#ifdef PRINTALL
	printf("Target = %f\n", target);
	printf("Output = %f\n", output);
	
	printf("Delta output sum = %f\n", delta_output_sum);
	#endif

	if (inputseq == 0) {
		inputseq0 = error;
	} else if (inputseq == 1) {
		inputseq1 = error;
	} else if (inputseq == 2) {
		inputseq2 = error;
	} else {
		inputseq3 = error;
	}
	averageerror = (inputseq0+inputseq1+inputseq2+inputseq3)/4;


	return delta_output_sum;

}

// Back calculate the output weights based ont he delta sum
int getNewOutputWeights(double *output_weight_matrix, double delta_output_sum, double *inner_layer_matrix) {
	for (int i=0; i<NUM_NODES; i++) {
		output_weight_matrix[i] = output_weight_matrix[i] - delta_output_sum * inner_layer_matrix[i] * (double)SPEED;
		#ifdef PRINTALL
		printf("New Output Weights = %f \n", output_weight_matrix[i]);
		#endif
	}
	return 0;
}

// Get the error of the inner sum
int getDeltaInnerSum(double delta_output_sum, double *output_weight_matrix, double *inner_layer_matrix, double *delta_inner_sum) {
	
	for (int i=0; i<NUM_NODES; i++) {
		delta_inner_sum[i] = delta_output_sum * output_weight_matrix[i] * applyActivationDerivative(inner_layer_matrix[i]);
		#ifdef PRINTALL
		printf("Delta Inner Sum = %f \n", delta_inner_sum[i]);
		#endif
	}

	return 0;
}

int getNewInnerWeights(double *input_weight_matrix, double *delta_inner_sum, double *inputs) {

	for (int i=0; i<NUM_NODES; i++) {
		for (int j=0; j<NUM_INPUTS; j++) {
			
			input_weight_matrix[i+j*NUM_NODES] = input_weight_matrix[i+j*NUM_NODES] - delta_inner_sum[i] * inputs[j] * (double)SPEED;  
			#ifdef PRINTALL
			printf("NEW Input Weight matrix =  %f \n", input_weight_matrix[i+j*NUM_NODES]);
			#endif

		}		
	}

	return 0;
}


// Offset inputs to avoid divide by zero
int offsetInputs(int *inputs, int num) {
}


void ClearScreen()
    {
    int n;
    for (n = 0; n < 10; n++)
      printf( "\n\n\n\n\n\n\n\n\n\n" );
    }

int main () {

	inputseq = 0;
	inputseq0 = 1;
	inputseq1 = 1;
	inputseq2 = 1;
	inputseq3 = 1;

	int num_input_weights = NUM_INPUTS * NUM_NODES;
	int num_output_weights = NUM_NODES;

	double inputs[NUM_INPUTS];
	double inner_layer_matrix[NUM_NODES];
	double delta_inner_sum[NUM_NODES];

	double input_weight_matrix[NUM_NODES * NUM_INPUTS];
	double output_weight_matrix[NUM_NODES];
	
	init_system();

	

	// Randomize weights for input weights and output weights matrix
	randomizeWeights(input_weight_matrix, NUM_NODES, NUM_INPUTS);
	randomizeWeights(output_weight_matrix, NUM_NODES, 1);

	for (int iter=0; iter<1000000; iter++)
	{
		//ClearScreen();
		getInputs(inputs, NUM_INPUTS);
		// Compute inner layer node values and output value
		computeInnerLayers(inputs, input_weight_matrix, inner_layer_matrix);
		double output = computeOutput(output_weight_matrix, inner_layer_matrix);

		double target = getTarget(inputs);

		// Error calculation and correction 
		double delta_output_sum = getDeltaOutputSum(target, output);
		getDeltaInnerSum(delta_output_sum, output_weight_matrix, inner_layer_matrix, delta_inner_sum);

		getNewOutputWeights(output_weight_matrix, delta_output_sum, inner_layer_matrix);
		getNewInnerWeights(input_weight_matrix, delta_inner_sum, inputs);
		

		if (iter%100 == 0) {
			printf("0 XOR 0    Target = 1.0   Predicted = %f\n", 1.0 - inputseq0);
			printf("0 XOR 1    Target = 0.0   Predicted = %f\n", inputseq1);
			printf("1 XOR 0    Target = 0.0   Predicted = %f\n", inputseq2);
			printf("1 XOR 1    Target = 1.0   Predicted = %f\n", 1.0 - inputseq3);
			printf("Average Error = %f\n", averageerror);
			printf("Epoch %i\n", iter);
			printf("\n");

			system("clear");
		}
	}
}