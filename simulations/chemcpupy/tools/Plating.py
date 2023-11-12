#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: Plating.py
@author: chrisarcadia
@description: methods to assist in plate (container) preparation
@created: Thu May 31 10:00:03 2018
"""

import numpy
import scipy.special 
from chemcpupy import Containers

# data related methods
    
def make_permutation(mydata,seed=12345):
    numpy.random.seed(seed)
    perm = numpy.random.permutation(mydata.size)
    antiperm = perm.argsort()
    permdata = mydata[perm]
    return perm,antiperm,permdata

def vector_to_blocks(vector, block_width):
    vector = numpy.array(vector);    
    blocks = [];
    padded = False;
    for indblock in range(0,len(vector),block_width):
        blocks.append(vector[indblock:indblock+block_width]);
    pad_length = numpy.ceil(vector.size/block_width)*block_width - vector.size;
    if pad_length>0: # if needed, add zero padding to make blocks equal in length
        blocks[-1] = numpy.append(blocks[-1],numpy.zeros(int(pad_length),dtype=blocks[-1].dtype));
        padded = True;
    return blocks, {'block_width':block_width, 'block_count':len(blocks), 'padded':padded, 'pad_length':pad_length, 'data_length':len(vector), 'bits_per_block':block_width}

def blocks_to_vector(blocks, vector_width):    
    vector = numpy.array(blocks).reshape(-1)
    return vector[0:vector_width]

def vector_to_subdivided_constant_weight_blocks(vector, block_width, block_weight): 
    # generate constant weight* blocks by subdividing the block (block_width) into divisions of length block_weight and placing a single bit as high in each subdivision
    # *weight is the number of 1's (high bits) in a block
    # note: due to the subdivison, this encoding does not span the entire space of constant-weight codes of length block_width
    vector = numpy.array(vector);
    blocks = [];
    padded = False;    
    valid_weight = (block_width % block_weight) == 0; # a weight is valid if it is an integer factor of the width (as block_weight is equal to the subdivision size)    
    if valid_weight:     
        states_per_subdivision = int(block_width / block_weight);
        bits_per_subdivision = int(numpy.ceil(numpy.log2(states_per_subdivision)));
        bits_per_block = bits_per_subdivision * block_weight;  
        pad_length = 0;
        # for each block
        for indblock in range(0,len(vector),bits_per_block):
            block = numpy.array([],dtype=vector.dtype); # initilize block
            vectorblock = vector[indblock:indblock+bits_per_block]; 
            pad_length = bits_per_block - len(vectorblock);
            if pad_length>0: # if needed, add zero padding to make blocks equal in length                
                vectorblock = numpy.append(vectorblock,numpy.zeros(int(pad_length),dtype=vectorblock.dtype));
                padded = True;            
            # for each subdivision (sub-block)   
            for inddiv in range(0,len(vectorblock),bits_per_subdivision):
                subblock = numpy.zeros(block_weight); # initilize sub-block                
                vectorsubblock = vectorblock[inddiv:(inddiv+bits_per_subdivision)]; 
                # convert to vectorsubpart binary value and use it as the index of the bit to toggle high in the sub-block (little endian)
                index = 0; 
                for b in range(0,bits_per_subdivision): 
                    index = index + vectorsubblock[b]*2**(bits_per_subdivision-b-1); # index is decimal representation of subblock vector read as little endian binary string
                subblock[index] = 1; 
                block = numpy.append(block,subblock);
            blocks.append(block);  
        return blocks, {'block_width':block_width, 'block_weight':block_weight, 'block_count':len(blocks),  'padded':padded, 'pad_length':pad_length, 'data_length':len(vector), 'bits_per_block':bits_per_block, 'subdivisions':block_weight}                      
    else:
        raise ValueError('Invalid weight given. Make sure block weight is an integer factor of block width.')                

def subdivided_constant_weight_blocks_to_vector(blocks, vector_width, block_weight): 
    blocks = numpy.array(blocks);
    vector = [];
    errors_found = 0; # tally of the number of found sub-block errors (where weight or block does not match block_weight)
    block_width = len(blocks[0]);    
    states_per_subdivision = int(block_width / block_weight);
    bits_per_subdivision = int(numpy.ceil(numpy.log2(states_per_subdivision)));
    for block in blocks:
        subblocks = numpy.reshape(block, (block_weight,states_per_subdivision));
        for subblock in subblocks:
            index = numpy.where(subblock==1)[0];
            if not len(index)==1:
                errors_found = errors_found + 1; #True;
            if index.size == 0:  # if no ones present select a one at random
                index = numpy.random.choice(range(0,states_per_subdivision));
            else: # else select at random of the ones found
                index = int(numpy.random.choice(index));
            bitsStr = numpy.binary_repr(index, width=bits_per_subdivision);
            bits = [int(x) for x in bitsStr];
            vector = vector + bits;    
    return vector[0:vector_width], errors_found;


def vector_to_codebook_constant_weight_blocks(vector, block_width, block_weight, book_size, seed=None): 
    # generate constant weight* blocks by using a pseudo-randomized codebook of size book_size (in bits) of constant weight block_weight and of word length block_width
    # *weight is the number of 1's (high bits) in a block and also the number of message bits (there are 2^block_weight states in the codebook)
    vector = numpy.array(vector);
    blocks = [];
    padded = False;    
    pad_length = 0;
    # generate codebook
    numpy.random.seed(seed);    
    num_codewords = numpy.power(2,book_size);
    bits_per_block = book_size; #block_weight;    
    max_num_codewords = scipy.special.comb(block_width,block_weight);     
    if num_codewords>max_num_codewords:
        raise ValueError('Requested number of codewords exceeds the maximum number of distinct codewords.')                
    codebook = ['']*num_codewords; #numpy.zeros((num_codewords,block_width),dtype=numpy.int8);
    wordbasis = ''.join(['0'] * (block_width - block_weight) + ['1'] * (block_weight)); 
    #min_codeword_distance = numpy.Inf; # min intra-codeword distance
    for n in range(0,num_codewords):      
        making_word = True;
        while making_word: # loop to ensure codewords are unique
            shuffling = numpy.random.permutation(block_width);
            codeword = ''.join([wordbasis[i] for i in shuffling]);
            if codeword not in codebook:
                making_word = False;
        codebook[n] = codeword;    
        #print('made ' + str(n) + 'th codeword')
    codebookArray = numpy.array([[int(x) for x in y] for y in codebook],dtype=numpy.int8);    
    # for each block
    for indblock in range(0,len(vector),bits_per_block):
        block = numpy.array([],dtype=vector.dtype); # initilize block
        vectorblock = vector[indblock:indblock+bits_per_block]; 
        pad_length = bits_per_block - len(vectorblock);
        if pad_length>0: # if needed, add zero padding to make blocks equal in length                
            vectorblock = numpy.append(vectorblock,numpy.zeros(int(pad_length),dtype=vectorblock.dtype));
            padded = True;            
        # for each subdivision (sub-block)       
        index = 0; 
        for b in range(0,bits_per_block): 
            index = index + vectorblock[b]*2**(bits_per_block-b-1); # index is decimal representation of subblock vector read as little endian binary string        
        block = codebookArray[index,:];
        blocks.append(block);  
    return blocks, {'block_width':block_width, 'block_weight':block_weight, 'book_size':book_size, 'block_count':len(blocks),  'padded':padded, 'pad_length':pad_length, 'data_length':len(vector), 'bits_per_block':bits_per_block, 'codebook':codebookArray}                      

def codebook_constant_weight_blocks_to_vector(blocks, vector_width, block_weight, codebook): 
    blocks = numpy.array(blocks);
    vector = [];
    block_distance = []; # distance of block to nearest codeword
    num_codewords = len(codebook);
    bits_per_block = int(numpy.floor(numpy.log2(num_codewords))); #block_weight;    
    for block in blocks:
        # compute Manhattan/Hamming distance to all codewords
        distance = numpy.inf + numpy.zeros(num_codewords);
        for n in range(1,num_codewords):
            distance[n] = numpy.sum(numpy.abs(codebook[n,:]-block)); # this is the Manhattan distance but for binary alphabet this is equivalent to Hamming distance          
        index = numpy.argmin(distance);  
        block_distance.append(distance[index]);          
        bitsStr = numpy.binary_repr(index, width=bits_per_block);
        bits = [int(x) for x in bitsStr];
        vector = vector + bits;    
    return vector[0:vector_width], block_distance;


def random_blocks(number_of_blocks, block_width, block_weight,seed=None):
    numpy.random.seed(seed);
    blocks = [];
    for n in range(0,number_of_blocks):
        new_block = numpy.array([0] * (block_width - block_weight) + [1] * (block_weight), dtype=numpy.int8);
        numpy.random.shuffle(new_block); # since we will not limit number_of_blocks by number of possible blocks, duplicates may occur
        blocks = blocks + [new_block];
    return blocks

# position list related methods

def update_free_positions(available,used):
    for pos in set(used):
        if pos in available:
            available.remove(pos);
    
def get_list_elements_by_indices(full_list,indices):
    return list(map(full_list.__getitem__, indices));

def sort_list_pair_by_first(X,Y):
    Xsorted = [];
    Ysorted = [];
    for (x,y) in sorted(zip(X,Y), key=lambda pair: pair[0]):
        Xsorted.append(x);
        Ysorted.append(y);
    return [Xsorted,Ysorted];

def sort_list_pair_by_first_then_second(X,Y):    
    # sort both by first list
    Xs,Ys = sort_list_pair_by_first(X,Y);       
    # sort like-valued bins of the first lit by their second list values 
    Xsorted = [];
    Ysorted = [];
    nbin = [];  
    xlast = not(Xs[0]);  
    N = len(Xs);
    Xs = Xs + [not(Xs[-1])];
    for n in range(0,N+1):        
        xcurr = Xs[n];
        if xcurr==xlast and n<N:
            nbin = nbin + [n];                
        else:
            xbin = get_list_elements_by_indices(Xs,nbin);
            ybin = get_list_elements_by_indices(Ys,nbin);
            ys,xs = sort_list_pair_by_first(ybin,xbin);    
            Xsorted = Xsorted + xs;
            Ysorted = Ysorted + ys;       
            nbin = [n];
        xlast = xcurr;        
    return [Xsorted,Ysorted];

def sort_by_source_then_destination(source,destination):
    return sort_list_pair_by_first_then_second(Containers.lettergrid_to_position(source), Containers.lettergrid_to_position(destination)); # for proper ordering we must ensure we are in position space and not letter space 

def sort_list_of_positions(positions): # assumed that positions are numerical tuple form and sorts by Y (first) and then X (second)
    return sorted(positions, key=lambda tup: (tup[0],tup[1]));
