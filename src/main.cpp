#include<stdio.h>
#include<iostream>
#include<math.h>

#include <unistd.h> 
#include <armadillo>

using namespace arma;

int main(int argc, char **argv)
{
	printf("Meshfree AD\n");

	/* initialize random seed*/
	srand (time(NULL));
	arma_rng::set_seed_random();

	

}