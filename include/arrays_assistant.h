/**
 * \file arrays_assistant.h
 * \author Alex Andriati
 * \date July 2021
 * \brief File with clash-protected basic auxiliar functions for arrays
 *
 * This file is supposed to be (absolutely) private for end users which
 * will consume the library. Moreover, the static keywords ensure clash
 * cannot happen with other libraries and the implementation cannot be
 * duplicated within this library. Only include if using all functions
 * otherwise warning will be prompt in compilation
 */

#ifndef ARRAYS_ASSISTANT_H
#define ARRAYS_ASSISTANT_H

#include "arrays.h"
#include <stdlib.h>
#include <stdio.h>


/** \brief Return fresh allocated real(double) array */
static Rarray
alloc_rarr(unsigned int array_size)
{
    Rarray ptr = (Rarray) malloc(array_size * sizeof(double));
    if (ptr == NULL)
    {
        printf("\n\nProblem in Rarray allocation\n\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}


/** \brief Return fresh allocated complex(double) array */
static Carray
alloc_carr(unsigned int array_size)
{
    Carray ptr = (Carray) malloc(array_size * sizeof(double complex));
    if (ptr == NULL)
    {
        printf("\n\nProblem in Carray allocation\n\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}


/** \brief Copy values from the first array to the second */
static void
carr_copy_values(unsigned int array_size, Carray from, Carray to)
{
    for (unsigned int i = 0; i < array_size; i++) to[i] = from[i];
}


/** \brief Copy values from the first array to the second */
static void
rarr_copy_values(unsigned int array_size, Rarray from, Rarray to)
{
    for (unsigned int i = 0; i < array_size; i++) to[i] = from[i];
}


#endif
