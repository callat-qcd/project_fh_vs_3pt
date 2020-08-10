//
//  main.c
//  energy_level_calculator
//
//  Created by Ivan Chernyshev on 3/22/20.
//  Copyright © 2020 Ivan Chernyshev. All rights reserved.
//
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Define a sort algorithm

int main(int argc, const char * argv[]) {
    // Take in number of pions and max pion momentum as command line arguments
    // Exit the program if the number of command-line arguments isn't 2
    printf("%i\n", argc);
    
    if(argc != 3)
    {
        printf("%s\n", "Error: energy_level_calculator must be called with 2 command-line arguments. The first is the number of pions, the second is the maximum magnitude of the momentum coordinate.");
        exit(EXIT_FAILURE);
    }
    
    int n_pi = atoi(argv[1]);
    int pmax = atoi(argv[2]);
    
    const int n_dim = 3; // Number of dimensions in a momentum array
    const double m_N = 0.668; // Nuclear mass in lattice units
    const double m_pi = 0.1885; // Pion mass in lattice units
    const double p_quantum = M_PI/12.; // The momentum-quantum in this lattice (of length 24)
    
    // Get the input for the number of pions and the maximum momentum coordinate
    char input[5];
    //printf("Insert number of pions (MUST BE AN INTEGER. ENTER 0 or a negative number to abort): ");
    //scanf("%s", input);
    //int n_pi = atoi(input); // Number of pions
    
    //printf("Insert maximal momentum index (MUST BE AN INTEGER. Enter a negative number to abort): ");
    //scanf("%s", input);
    //int pmax = atoi(input); // Maximum momentum coordinates
    
    printf("Number of pions: %i; maximum momentum coordinate: %i\n", n_pi, pmax);
    
    // Abort the program if a non-positive npi or a negative pmax was entered
    if(n_pi <= 0 || pmax < 0)
    {
        printf("Exiting.\n");
        exit(EXIT_FAILURE);
    }
    
    // Warning if an attempt to make over a billion data samples is made
    if(pow(2*pmax+1, 3*(n_pi -1))*(pmax + 1)*(pmax + 2)*(pmax + 3)/6. > 1000000000)
    {
        printf("WARNING: ATTEMPT TO MAKE OVER 1 BILLION DATA SAMPLES. YOUR DATA SIZE COULD EASILY REACH TENS OF GIGABYTES OR ORDERS OF MAGNITUDE MORE. IF YOU PROCEED, THE CREATOR IS NOT RESPONSIBLE FOR ANY SYSTEM PROBLEMS. PROCEED? (y/n)");
        scanf("%s", input);
        if(strcmp(input, "y") != 0)
        {
            printf("Exiting.\n");
            exit(EXIT_FAILURE);
        }
        
    }
        
    
    // Failsafe: crashes the program if the total number of choices exceeds the maximum value representable by long long int (9,223,372,036,854,775,807)
    if(pow(2*pmax+1, 3*(n_pi -1))*(pmax + 1)*(pmax + 2)*(pmax + 3)/6. > LLONG_MAX)
    {
        printf("Congratulations. You just tried to create more data than a signed long long could count with just positive numbers (over 9 exabytes). Happy now?\nEnding program.");
        exit(EXIT_FAILURE);
    }
    
    // Define quantities
    // Storage of results
    unsigned long long int num_nucleon_p = (pmax + 1)*(pmax + 2)*(pmax + 3)/6; // Number of possible nucleon momenta, calculated analyticaly
    unsigned long long int n_pchoices = 1; // Number of possible pion momentum combinations, calculated analytically
    for(int i = 0; i < 3*(n_pi - 1); i++)
        n_pchoices *= 2*pmax+1;
    unsigned long long int n_choices = num_nucleon_p*n_pchoices; // Total number of possible momentum picks, calculated analytically
    printf("Number of nucleon momentum  choices: %llu\n", num_nucleon_p);
    printf("Number of pion momentum  choices: %llu\n", n_pchoices);
    printf("Number of total momentum choices: %llu\n", n_choices);
    
    //printf("Number of nucleon momentum  choices: %12.2f\n", (pmax + 1.)*(pmax + 2.)*(pmax + 3.)/6.);
    //printf("Number of pion momentum  choices: %12.2f\n", pow(2*pmax+1, 3*(n_pi -1)));
    //printf("Number of total momentum choices: %12.2f\n", pow(2*pmax+1, 3*(n_pi -1))*(pmax + 1)*(pmax + 2)*(pmax + 3)/6.);
    
    FILE *fptr_nuc; // pointer to the file being written to
    char filename[35];
    sprintf(filename, "energy_levels_npi_%i_pmax_%i.txt", n_pi, pmax);
    fptr_nuc = fopen(filename, "w");
    
    double E; // The holder for the energy
    
    // Print out a header for the file
    fprintf(fptr_nuc, "Energy|p_Nx, p_Ny, p_Nz|");
    for(int i = 0; i < n_pi-1; i++)
        fprintf(fptr_nuc, "p_pi%ix, p_pi%iy, p_pi%iz; ", i, i, i);
    
    fprintf(fptr_nuc, "p_pi%ix, p_pi%iy, p_pi%iz\n ", n_pi, n_pi, n_pi);
    
    // Record the default lengths of each line in the file, allowing for up to 2 digits in the ints in the neader, up to 3 for the ints in the main data stash, and up to 4 for the final pion
    
    
    // Loop over nuclear momentum indices
    int choice_index = 0; // The current entry-index for the choices of p_N, p_pi, and E to be inserted
    
    for (int nx = 0; nx <= pmax; nx++)
        for(int ny = 0; ny <= nx; ny++)
            for(int nz = 0; nz <= ny; nz++)
            {
                // Calculate energy contribution from the nucleon
                E = sqrt(m_N*m_N + (p_quantum*nx)*(p_quantum*nx) + (p_quantum*ny)*(p_quantum*ny) + (p_quantum*nz)*(p_quantum*nz));
                for(int i = 0; i < n_pi - 1; i++)
                    E += sqrt(m_pi*m_pi + (p_quantum*(-pmax))*(p_quantum*(-pmax)) + (p_quantum*(-pmax))*(p_quantum*(-pmax)) + (p_quantum*(-pmax))*(p_quantum*(-pmax)));
                
                E += sqrt(m_pi*m_pi + p_quantum*((n_pi - 1)*pmax-nx)*(p_quantum*((n_pi - 1)*pmax-nx)) + p_quantum*((n_pi - 1)*pmax-ny)*(p_quantum*((n_pi - 1)*pmax-ny)) + p_quantum*((n_pi - 1)*pmax-nz)*(p_quantum*((n_pi - 1)*pmax-nz)));
                
                fprintf(fptr_nuc, "%4.4f|%i,%i,%i|", E, nx, ny, nz);
                for(int i = 0; i < n_pi - 1; i++)
                    fprintf(fptr_nuc, "%i,%i,%i;", -pmax, -pmax, -pmax);
                
                fprintf(fptr_nuc, "%i,%i,%i\n", (n_pi - 1)*pmax-nx, (n_pi - 1)*pmax-ny, (n_pi - 1)*pmax-nz);
                
            }
    fclose(fptr_nuc);
    printf("Writing of nuclear momenta and lowest pion momenta (all coordinates %i) complete.\n",-pmax);
    
    // This is the pion momentum section
    
    
    FILE *fptr_pion_r; // pointer to the file being written to, meant to read the file
    FILE *fptr_pion_a; // pointer to the file being written to, meant to read append to the file
    
    char header_read[26*n_pi+27]; // Array for reading in the header-line
    char lineread[12*n_pi+23]; // Array for reading in a line in the main body of the data
    char* line_split; // Array for reading in a section of the line read, separated by the '|' character
    char* mom_split; // Array for reading in a pion momentum
    char momentum_coordinate[12] = "%i,%i,%i"; // Template for creating a momentum-coordinate
    char newline_template[12*n_pi+23]; // The template of the new line to be written to the new file
    double E_pi_standard = sqrt(m_pi*m_pi + 3*(p_quantum*(-pmax))*(p_quantum*(-pmax))); //This variable is for calculating the adjustment to the energy of the system from changing a pion momentum
    //double E_pi_adj;
    //int n_pix;// This and the next 3 coord
    //int n_piy;
    //int n_piz;
    char* last_pion_momentum; // This char pointer is there to point to the location of the last pion momentum
    char* nuc_momentum; // This char pointer points to the
    int n_lastreadx; // this and the next 3 coordinates are for the last pion momentum in the line being read
    int n_lastready;
    int n_lastreadz;
    int n_lastwrittenx; // this and the next 3 coordinates are for the last pion momentum in the line being written
    int n_lastwritteny;
    int n_lastwrittenz;
    long long int num_of_lines_read = num_nucleon_p;
    for(int i = 0; i < n_pi - 1;i++) // Iterate over the momenta of each pion except the last
    {
        fptr_pion_r = fopen(filename, "r");
        fptr_pion_a = fopen(filename, "a");
        
        fgets(header_read, 26*n_pi+26, fptr_pion_r); // use fgets to read off the header
        
        // Read in (2*pmax + 1)^i lines (pmax + 1)(pmax + 2)(pmax + 3)/6
        for(long long int j = 0; j < num_of_lines_read; j++)
        {
            fgets(lineread, 12*n_pi+22, fptr_pion_r); // read in line
            
            // Get the momentum of the last pion
            int k = strlen(lineread);
            for(; k >= 0; k--)
            {
                //printf("%c\n", lineread[k]);
                if(lineread[k] == ';')
                    break;
            }
            
            if(k < 0) // If the 'get the momentum of the last pion' subroutine failed, that means that there are no pion momenta in the current data set. In this case, the current data set is skipped
            {
                printf("Error: unable to find pion momenta in the datastash. Skipping current loop iteration\n");
                continue;
            }
            
            last_pion_momentum = lineread + k + 1; // Assign the last pion momentum pointer to the location right after the last semicolon
            mom_split = strtok(last_pion_momentum, ",");// Load data of the last pion momentum into n_lastreadx, n_lastready, n_lastreadz
            n_lastreadx = atoi(mom_split);
            mom_split = strtok(NULL, ",");
            n_lastready = atoi(mom_split);
            mom_split = strtok(NULL, ",");
            n_lastreadz = atoi(mom_split);
            // Since strtok replaces all of its delimiters with null chars, I need to put the delimiters back in so that the strings can be reused
            last_pion_momentum[strlen(last_pion_momentum)] = ',';
            last_pion_momentum[strlen(last_pion_momentum)] = ',';
            
            // Get the energy
            line_split = strtok(lineread, "|");
            E = atof(line_split);
            line_split[strlen(line_split)] = '|';
            
            // Loop over all of the possible coordinates of the relevant pions
            for(int n_pix = -pmax; n_pix <= pmax; n_pix++)
                for(int n_piy = -pmax; n_piy <= pmax; n_piy++)
                    for(int n_piz = -pmax; n_piz <= pmax; n_piz++)
                    {
                        if(n_piz == -pmax && n_piy == -pmax && n_pix == -pmax)
                            continue;
                        // Calculate the correction to the last pion's momentum
                        n_lastwrittenx = n_lastreadx - pmax - n_pix;
                        n_lastwritteny = n_lastready - pmax - n_piy;
                        n_lastwrittenz = n_lastreadz - pmax - n_piz;
                        
                        // Write the corrected energy to the file
                        fprintf(fptr_pion_a, "%4.4f|", E - E_pi_standard + sqrt(m_pi*m_pi + (n_pix*p_quantum)*(n_pix*p_quantum) + (n_piy*p_quantum)*(n_piy*p_quantum) + (n_piz*p_quantum)*(n_piz*p_quantum)) - sqrt(m_pi*m_pi + (n_lastreadx*p_quantum)*(n_lastreadx*p_quantum) + (n_lastready*p_quantum)*(n_lastready*p_quantum) + (n_lastreadz*p_quantum)*(n_lastreadz*p_quantum)) + sqrt(m_pi*m_pi + (n_lastwrittenx*p_quantum)*(n_lastwrittenx*p_quantum) + (n_lastwritteny*p_quantum)*(n_lastwritteny*p_quantum) + (n_lastwrittenz*p_quantum)*(n_lastwrittenz*p_quantum)));
                        
                        // Copy over the nucleon momentum in the written file
                        line_split = strtok(lineread, "|");
                        line_split = strtok(NULL, "|");
                        fprintf(fptr_pion_a, "%s|", line_split);
                        line_split = strtok(NULL, "|");
                        mom_split = strtok(line_split, ";");
                    
                        // Now, copy over the pion momenta that are not to be replaced
                        for(int l = 0; l < i; l++)
                        {
                            fprintf(fptr_pion_a, "%s;", mom_split);
                            mom_split = strtok(NULL, ";");
                        }
                        // Copy over the re-written pion momentum into the file
                        fprintf(fptr_pion_a, "%i,%i,%i;", n_pix,n_piy,n_piz);
                        mom_split = strtok(NULL, ";");
                        
                        for(int l = i + 1; l < n_pi - 1; l++)
                        {
                            fprintf(fptr_pion_a, "%s;", mom_split);
                            mom_split = strtok(NULL, ";");
                        }
                        
                        // And finally, the last pion's momentum
                        fprintf(fptr_pion_a, "%i,%i,%i\n", n_lastwrittenx, n_lastwritteny, n_lastwrittenz);
                        
                        // Finally, reset the delimiters that strtok set to NULL back to what they were before
                        lineread[strlen(lineread)] = '|';
                        lineread[strlen(lineread)] = '|';
                        for(int l = 0; l < n_pi - 1; l++)
                            lineread[strlen(lineread)] = ';';
                    }
            
            // Begin constructing the template for the new line
            //strcat(newline_template, "%4.4f|");
            
            // Load the nucleon momentum into the template
            //line_split = strtok(NULL, "|");
            //strcat(newline_template, line_split);
            //strcat(newline_template, "|");
            
            // Start on the pion momenta. Simply load-in all the pion momenta excpet the ith one and the last one (which will need re-adjustment
            //line_split = strtok(NULL, "|");
            //printf("Secondtry: %s", line_split);
            //line_split = strtok(NULL, "|");
            //printf("Thirdtry: %s", line_split);
            //mom_split = strtok(line_split, ";");
            //printf(" Momvector_1:(%s)", mom_split);
            
            //mom_split = strtok(NULL, ";");
            //printf(" Momvector_2:(%s)", mom_split);
            //mom_split = strtok(NULL, ";");
            //printf(" Momvector_3:(%s)", mom_split);
            //mom_split = strtok(line_split, ";");
            //printf(" Momvector_ret:(%s)", mom_split);
            
            // Now, test out the semicolon-delimiter strtok
            //for(int i = 0; i < n_pi; i++)
                //printf("Semicolon-delimiter: %s", strtok(line_split, ";"));
            
            //line_split = strtok(lineread, "|");
            //printf("Retry: %s", line_split);
            
            //printf("\n");
            
            
            //memset(newline_template,0,strlen(newline_template)); // Reset the memory of the newline_template string
        }
        
        
        fclose(fptr_pion_r);
        fclose(fptr_pion_a);
        
        num_of_lines_read *= (2*pmax + 1);
        num_of_lines_read *= (2*pmax + 1);
        num_of_lines_read *= (2*pmax + 1);// Recalculate the number of lines to read for the next iteration
    }
    
    // Arrays containing information about the choices possible
    //int p_N[n_choices][n_dim];
    //int p_pi[n_choices][n_pi][n_dim];
    //double E[n_choices];
    //for(int i = 0; i < n_choices; i++)
        //E[i] = 0;
    
    //int Nuc[n_dim]; // The nuclear momentum vector
    //int p_cumulative[n_dim]; // The cumulative momentum vector
    //int p_pis[n_dim*(n_pi-1)]; // The pion momentum vectors, except one, which is set by the "sum momentum zero" rule
    
// Loop over nuclear momentum indices
/*
    int choice_index = 0; // The current entry-index for the choices of p_N, p_pi, and E to be inserted
    
    for (int nx = 0; nx <= pmax; nx++)
        for(int ny = 0; ny <= nx; ny++)
            for(int nz = 0; nz <= ny; nz++)
            {
                p_N[choice_index][0] = nx;
                p_N[choice_index][1] = ny;
                p_N[choice_index][2] = nz;
                
                // Set all pion momenta numbers other than that of the last pion to the 'default' value (i.e. -pmax)
                // The last pion is the additive inverse of the sum of the all the other momenta
                // Simulataneously use the regular relativistic energy formula to calculate the total energy of the system
                p_pi[choice_index][n_pi - 1][0] = -p_N[choice_index][0];
                p_pi[choice_index][n_pi - 1][1] = -p_N[choice_index][1];
                p_pi[choice_index][n_pi - 1][2] = -p_N[choice_index][2];
                E[choice_index] = sqrt(m_N*m_N + (p_quantum*p_N[choice_index][0])*(p_quantum*p_N[choice_index][0]) + (p_quantum*p_N[choice_index][1])*(p_quantum*p_N[choice_index][1]) + (p_quantum*p_N[choice_index][2])*(p_quantum*p_N[choice_index][2]));
                
                for(int j = 0; j < n_pi - 1; j++)
                {
                    p_pi[choice_index][j][0] = -pmax;
                    p_pi[choice_index][j][1] = -pmax;
                    p_pi[choice_index][j][2] = -pmax;
                    
                    p_pi[choice_index][n_pi - 1][0] -= p_pi[choice_index][j][0];
                    p_pi[choice_index][n_pi - 1][1] -= p_pi[choice_index][j][1];
                    p_pi[choice_index][n_pi - 1][2] -= p_pi[choice_index][j][2];
                    
                    E[choice_index] += sqrt(m_pi*m_pi + (p_quantum*p_pi[choice_index][j][0])*(p_quantum*p_pi[choice_index][j][0]) + (p_quantum*p_pi[choice_index][j][1])*(p_quantum*p_pi[choice_index][j][1]) + (p_quantum*p_pi[choice_index][j][2])*(p_quantum*p_pi[choice_index][j][2]));
                }
                
                E[choice_index] += sqrt(m_pi*m_pi + (p_quantum*p_pi[choice_index][n_pi - 1][0])*(p_quantum*p_pi[choice_index][n_pi - 1][0]) + (p_quantum*p_pi[choice_index][n_pi - 1][1])*(p_quantum*p_pi[choice_index][n_pi - 1][1]) + (p_quantum*p_pi[choice_index][n_pi - 1][2])*(p_quantum*p_pi[choice_index][n_pi - 1][2]));
                
                choice_index++;
            }
    
    // For each of the other momentum indices of the first npi - 1 pions, loop through all momentum indices between -pmax (not inclusive) and pmax (inclusive), with all other elements being cloned from the previous (1/6)(pmax + 1)(2*pmax + 1)^i where i is the index of the current momentum coordinate
    // The nucleon momentum are simply cloned, as described above
    for(int i = 0; i < 3*(n_pi - 1); i++)
    {
        for(int pj = -pmax + 1; pj <= pmax; pj++)
            for(int j = 0; j < num_nucleon_p*pow(2*pmax + 1, i); j++)
            {
                // Nucleon momentum
                p_N[choice_index][0] = p_N[j][0];
                p_N[choice_index][1] = p_N[j][1];
                p_N[choice_index][2] = p_N[j][2];
                
                
                // Pion momenta
                // Simultaneously calculate the last pion's momentum
                p_pi[choice_index][n_pi-1][0] = -p_N[choice_index][0];
                p_pi[choice_index][n_pi-1][1] = -p_N[choice_index][1];
                p_pi[choice_index][n_pi-1][2] = -p_N[choice_index][2];
                
                for(int k = 0; k < i; k++)
                {
                    p_pi[choice_index][k/3][k % 3] = p_pi[j][k/3][k % 3];
                    p_pi[choice_index][n_pi-1][k % 3] -= p_pi[choice_index][k/3][k % 3];
                }
                p_pi[choice_index][i/3][i % 3] = pj;
                p_pi[choice_index][n_pi-1][i % 3] -= pj;
                for(int k = i + 1; k < 3*(n_pi - 1);k++)
                {
                    p_pi[choice_index][k/3][k % 3] = p_pi[j][k/3][k % 3];
                    p_pi[choice_index][n_pi-1][k % 3] -= p_pi[choice_index][k/3][k % 3];
                }
                
                // Now, calculate the energy of this choice
                E[choice_index] = sqrt(m_N*m_N + (p_quantum*p_N[choice_index][0])*(p_quantum*p_N[choice_index][0]) + (p_quantum*p_N[choice_index][1])*(p_quantum*p_N[choice_index][1]) + (p_quantum*p_N[choice_index][0])*(p_quantum*p_N[choice_index][2]));
                for(int k = 0; k < n_pi - 1; k++)
                    E[choice_index] += sqrt(m_pi*m_pi + (p_quantum*p_pi[choice_index][k][0])*(p_quantum*p_pi[choice_index][k][0]) + (p_quantum*p_pi[choice_index][k][1])*(p_quantum*p_pi[choice_index][k][1]) + (p_quantum*p_pi[choice_index][k][2])*(p_quantum*p_pi[choice_index][k][2]));
                choice_index++;
            }
    }
    
    
    for (int nx = 0; nx <= pmax; nx++)
        for(int ny = 0; ny <= nx; ny++)
            for(int nz = 0; nz <= ny; nz++)
            {
                Nuc[0] = nx;
                Nuc[1] = ny;
                Nuc[2] = nz;
                for (int i = 0; i < n_pchoices; i++)
                {
                    int j = 0;
                    int current_n = i;
                    while(current_n != 0) // Treats the possibilites of pion momenta as an "all values of a (2*pmax+1)-base 3(n_pi - 1)-digit number" problem
                    {
                        //printf("Value of current digit: %5.2f; inted: %i", pow((2*pmax+1), 3*(n_pi-1)-1-j), (int) pow((2*pmax+1), 3*(n_pi-1)-1-j));
                        int current_digitval = (int) pow((2*pmax+1), 3*(n_pi-1)-1-j);// The value of a unit in the current "digit" being filled
                        printf("Current pass: Nuclear momenum: (%i, %i, %i); P-choice index: %i; Current n: %i; Current digit: %i; order of current digit: %i\n", Nuc[0], Nuc[1], Nuc[2], i, current_n, current_digitval, j);
                        //printf("Next digit assignment: %2.2d \n\n", current_n/current_digitval - pmax);
                        p_pis[j] = current_n/current_digitval - pmax; //Fill the current digit, decrease the digit by pmax to cover both positive and negative momentum coordinates
                        
                        
                        // Increment up to the next digit
                        current_n = current_n % current_digitval;
                        j++;
                    }
                    // After all the "digits" required to be nonzero to represent the "number" are found, set all remaining "digits" to zero (i.e. minimum momentum value)
                    for(;j<3*(n_pi-1);j++)
                        p_pis[j] = -pmax;
                    
                    // Now, find the value of p_cumulative
                    // Simultaneously, fill out the choice arrays
                    for(int k = 0; k < n_dim; k++)
                    {
                        p_cumulative[k] = Nuc[k];
                        p_N[choice_index][k] = Nuc[k];
                    }
                    E[choice_index] = sqrt(m_N*m_N + (nx*p_quantum)*(nx*p_quantum) + (ny*p_quantum)*(ny*p_quantum) + (nz*p_quantum)*(nz*p_quantum));
                   
                    // After adding in the nuclear momentum, add in all the pion momenta
                    for(int k = 0; k < n_pi-1; k++)
                    {
                        //printf("Max p_pis size: %i; attempted access indices: %i, %i, %i; values: %i, %i, %i", n_dim*(n_pi-1), n_dim*k, n_dim*k+1, n_dim*k+2, p_pis[n_dim*k], p_pis[n_dim*k+1], p_pis[n_dim*k+2]);
                        p_cumulative[0] += p_pis[n_dim*k];
                        p_cumulative[1] += p_pis[n_dim*k+1];
                        p_cumulative[2] += p_pis[n_dim*k+2];
                        E[choice_index] += sqrt(m_pi*m_pi + (p_pis[n_dim*k]*p_quantum)*(p_pis[n_dim*k]*p_quantum) + (p_pis[n_dim*k+1]*p_quantum)*(p_pis[n_dim*k+1]*p_quantum) + (p_pis[n_dim*k+2]*p_quantum)*(p_pis[n_dim*k+2]*p_quantum));
                        
                        p_pi[choice_index][k][0] = p_pis[n_dim*k];
                        p_pi[choice_index][k][1] = p_pis[n_dim*k + 1];
                        p_pi[choice_index][k][2] = p_pis[n_dim*k + 2];
                    }
                    
                    // The final part of p_pi and is the additive inverse of the cumulative momentum
                    p_pi[choice_index][n_pi-1][0] = -p_cumulative[0];
                    p_pi[choice_index][n_pi-1][1] = -p_cumulative[1];
                    p_pi[choice_index][n_pi-1][2] = -p_cumulative[2];
                    E[choice_index] += sqrt(m_pi*m_pi + (p_cumulative[0]*p_quantum)*(p_cumulative[0]*p_quantum) + (p_cumulative[1]*p_quantum)*(p_cumulative[1]*p_quantum) + (p_cumulative[2]*p_quantum)*(p_cumulative[2]*p_quantum));
                    
                    printf("Current entry: Choice index: %i Energy: %1.4f Nucleon momentum: (%i, %i, %i) Pion momenta: ", choice_index, E[choice_index], p_N[choice_index][0], p_N[choice_index][1], p_N[choice_index][2]);
                    for(int l = 0; l < n_pi; l++)
                        printf("(%i, %i, %i), ", p_pi[choice_index][l][0], p_pi[choice_index][l][1], p_pi[choice_index][l][2]);
                    
                    printf("\n\n");
                    choice_index++;
                    
                }
                
                
            }
    */
    // Now, use merge-sort algorithm to sort the energies
    
    // Now, print out the possible options for the excited states
    /*
    for(int i = 0; i < n_choices; i++)
    {
        printf("\nEnergy of state: %2.3f; Momenta: nucleon: (%i, %i, %i); pions: ", E[i], p_N[i][0], p_N[i][1], p_N[i][2]);
        
        for(int j = 0; j < n_pi; j++)
            printf("(%i, %i, %i),", p_pi[i][j][0], p_pi[i][j][1], p_pi[i][j][2]);
        
    }
    printf("\n");
    */
    printf("Done producing all combinations.\n");
    return 0;
}
