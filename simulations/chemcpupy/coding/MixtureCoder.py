#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:04:32 2018

@author: jacobrosenstein
"""


import numpy as np
import chemcpupy as ccpu


class MixtureCoder:
    
    def __init__(self,num_molecules=32,num_mixtures=8,block_size=(4,8,4)):
        """ A MixtureCoder performs random block coding.
        num_molecules:   the number of molecules available in the library
        num_mixtures:    the number of mixtures available
        block_size:  3 dimensions:   (input bits, output molecules, output mixtures)
        e.g. (4,4,8) will code every 4 bits into a 4x8=32 bit block
        """
        assert(np.mod(num_molecules,block_size[0])==0)
        assert(np.mod(num_molecules,block_size[1])==0)
        assert(np.mod(num_mixtures,block_size[2])==0)
        assert(block_size[0]<=8)
        
        self._num_molecules=num_molecules
        self._num_mixtures=num_mixtures
        self._block_size=block_size

        np.random.seed(100)
        self._code_book = []
        for i in range(2**self._block_size[0]):
            self._code_book.append(np.random.randint(0,2,self._block_size[1:]))

        #self.load_data(np.empty(0,np.uint8))        

        
    def load_data(self,data):
        self._data = data
        self._data_uint8 = data.view(np.uint8)
        self._num_bytes = len(self._data_uint8)
        
        self._data_binary = np.unpackbits(self._data_uint8)
        self._num_bits = len(self._data_binary)
        
        coded_bits = self._num_molecules*self._num_mixtures*self._block_size[0]/np.prod(self._block_size[1:])
        print("%d bits loaded, %d will be encoded.\n" % (self._num_bits,coded_bits))
        if self._num_bits < coded_bits:
            raise Exception('not enough data bits supplied')
        elif self._num_bits > coded_bits:
            print("\nWARNING:  too much data supplied. Only the start of the message will be encoded.\n")
        
        self._coded_binary = np.zeros( (self._num_molecules,self._num_mixtures), dtype=np.uint8 )
        blockindices=np.arange(self._block_size[0],dtype=np.int)
        for i in range(0,self._num_molecules,self._block_size[1]):
            for j in range(0,self._num_mixtures,self._block_size[2]):
                codeindex = np.right_shift(np.packbits(self._data_binary[blockindices])[0],
                                           8-self._block_size[0])

                if 0:
                    print(self._data_binary[blockindices][0],'   ',codeindex)
                    for b in self._data_binary[blockindices][1:]:
                        print(b)
                
                blockoutput = self._code_book[codeindex]
                self._coded_binary[ i:(i+self._block_size[1]), j:(j+self._block_size[2]) ] = blockoutput
                blockindices += self._block_size[0]


    def decode_block(self,block):
        """ Decodes the block into the original data.
        """
        decoded_data = np.zeros_like(self._data_binary)
        self._decoded_distances = np.zeros( shape=(len(self._data_binary),len(self._code_book)) )
        
        blockindices=np.arange(self._block_size[0],dtype=np.int)
        for i in range(0,self._num_molecules,self._block_size[1]):
            for j in range(0,self._num_mixtures,self._block_size[2]):
                thisblock = block[ i:(i+self._block_size[1]), j:(j+self._block_size[2]) ]
                dist = np.zeros(len(self._code_book))
                for c in range(len(self._code_book)):
                    dist[c] = np.linalg.norm(thisblock - self._code_book[c],ord=2)

                #print('min distance',min(dist),np.argmin(dist),np.unpackbits(np.uint8(np.argmin(dist)))[-self._block_size[0]:])
                self._decoded_distances[blockindices,:] = dist
                decoded_data[blockindices] = np.unpackbits(np.uint8(np.argmin(dist)))[-self._block_size[0]:]
                
                blockindices += self._block_size[0]
        print('decoded shape',decoded_data.shape)
        return decoded_data

        
    def data_to_mixtures(self,mylibrary):
        """ Produces a MixtureList representing the data, from a provided CompoundList library.
        """        
        assert(len(mylibrary)>=self._num_molecules)        
        
        data_mixtures = ccpu.MixtureList()
        for m in range(self._num_mixtures):
            cl = mylibrary.__class__()
            cl.add_compounds(compound_list=[x for index,x in enumerate(mylibrary) if (index<self._num_molecules and self._coded_binary[index,m])])
            data_mixtures.add_mixture('data_mix_%d'%(m,), cl)
        return data_mixtures 


    def codebook_to_mixtures(self,mylibrary):
        """ Produces a MixtureList representing the possible codebook combinations, 
        from a provided CompoundList library. Note that the codebook may have several possible alignments
        with the molecules, and these are all included.
        """                
        codebook_mixtures = ccpu.MixtureList()
        for cb in range(len(self._code_book)):
            for k in range(int(self._num_molecules/len(self._code_book[0]))):
                cl = mylibrary.__class__()
                cl.add_compounds(compound_list=[x for index,x in enumerate(mylibrary) if int(index/len(self._code_book[0]))==k and self._code_book[cb][np.mod(index,len(self._code_book[0]))]])
                codebook_mixtures.add_mixture('codebook_mix_%d_%d'%(cb,k,), cl)
        return codebook_mixtures


    def mixtures_to_data(self,mylibrary,mymixtures):
        """ Recovers the binary data from a MixtureList representing the data, by comparing against 
        a CompoundList library.
        """
        
        m_index=0
        recovered_binary_data = np.zeros_like(self._coded_binary)
        for m in mymixtures:
            c_index=0
            if m_index<self._num_mixtures:
                for c in mylibrary:
                    if c_index<self._num_molecules:
                        if [c,] in m['compound_list']:
                            recovered_binary_data[c_index,m_index]=1
                        c_index+=1
                m_index+=1

        return recovered_binary_data
        
        

if __name__=='__main__':
    
    print('\n\n')
    
    np.random.seed(4)
    mydata = np.random.randint(0,100,4,dtype=np.uint8)
    
    x = MixtureCoder(block_size=(2,8,4))
    x.load_data(mydata)
    
    print('\n\n')
    print('code book')
    for i in range(len(x._code_book)):
        print('\n')
        print(i,' ',x._code_book[i][0])
        for c in x._code_book[i][1:]:
            print('   ',c)
    
    
    print('\n')
    print('data',x._data)
    print('\n')

    print('data             coded')
    for b,c in zip(x._data_binary,x._coded_binary):
        print(' ',b,'         ',c)







