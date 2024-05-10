# -*- coding: utf-8 -*- 
'''
Python 3.10 
Created Apr 30 2023 
Jasper Bellefeuille - belle172@umn.edu 
Repository: Human_Genome_Analysis/file_processing.py 

This script contains functions for doing text processing of files 
''' 

# Return the data from a tab seperated file as a matrix 
def get_tsv_matrix(filename): 
    file = open( filename, encoding = "utf8" ) 
    matrix = [] 
    for line in file: 
        if line[0] != '#': # go past comment lines 
            matrix.append(line.strip().split('\t')) 
    file.close() 
    return matrix 

# change matrix to only contain the specified columns 
def slim_matrix(matrix, columns): 
    slim_m = [] 
    for row in matrix: 
        slim_row = [] 
        for i in columns: 
            slim_row.append(row[i]) 
        slim_m.append(slim_row) 
    return slim_m 

# takes in a matrix and the key column containing the int or rsid to sort the matrix 
# returns the matrix with the specified column put as the first value in each row 
def rsid_matrix(matrix, index): 
    column = [] 
    for row in matrix: 
        rsid = row.pop(index) 
        if 'rs' in rsid: # accept 'rs<int>' or '<int>' as rsid 
            rsid = rsid[2:] 
        try: 
            output = [int(rsid)] + row 
            column.append( output ) 
        except ValueError: 
            1 
    return column 

