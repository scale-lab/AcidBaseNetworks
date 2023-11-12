#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 13:11:52 2018

@author: jacobrosenstein
"""

from lxml import etree
from time import sleep
import numpy as np
import csv

import chemcpupy.tools.Containers as Containers
#from chemcpupy.tools.misc import *

# ============================================================

class Task:
    
    def __init__(self,*args,**kwargs):
        self._task_type = None
        self._required_params = None
        self._params = None
        
    def _check_params(self):
        if not set(self._required_params).issubset(self._params):
            raise Exception('some required parameters are not present')
            
    def run(self,**kwargs):
        raise Exception('cannot run an undefined task')

    def __str__(self):
        return self._task_type

    def has_param(self, param):
        return param in self._params
        
    def get_param(self, param):
        return self._params.get(param)


# ============================================================

class PauseTask(Task):
    
    def __init__(self,*args,**kwargs):
        self._task_type = 'pause'        
        self._required_params = ('pause_time')
        
        self._params = kwargs
        self._check_params()

    def run(self,verbose=True):
        sleep(self._params['pause_time'])
        
    
# ============================================================

class TransferTask(Task):
    
    def __init__(self,*args,**kwargs):
        """ A TransferTask performs a list of liquid transfers 
        from one WellPlate to another.
        
        Kwargs:
            from_plate (WellPlate)
            from_positions (list of tuples)
            to_plate (WellPlate)
            to_positions (list of tuples)
            transfer_volumes (list of floats)
            pre_transfer_delays (optional)            
            transfer_group_label (optional)
            repeat_each   (optional)
        
        """
        self._task_type = 'transfer'        
        self._required_params = ('from_plate',
                                 'from_positions',
                                 'to_plate',
                                 'to_positions',
                                 'transfer_volumes')                    

        self._params = kwargs
        self._check_params()
        
    
    def unpack_transfers(self,rotate_dest_180=False):
        from_plate = self._params['from_plate']
        from_positions = self._params['from_positions']
        to_plate = self._params['to_plate']
        to_positions = self._params['to_positions']
        transfer_volumes = self._params['transfer_volumes']
        repeat_each = self._params.get('repeat_each',1)       # optional
        repeat_all = self._params.get('repeat_all',1)       # optional
        pre_transfer_delays = self._params.get('pre_transfer_delays',[0]) # optional
        
        # optional, runs assuming the destination plate is rotated 180 degrees
        if rotate_dest_180:
            print('SIMULATING 180 DEGREE ROTATED DESTINATION')
            to_positions = Containers.rotate_positions_180(to_positions,to_plate)
        
        # if position/volume is given as a single item, change to a single-element list
        if type(from_positions) in (tuple,str):
            from_positions = [from_positions,]
        if type(to_positions) in (tuple,str):
            to_positions = [to_positions,]
        if type(transfer_volumes) is int:
            transfer_volumes = [transfer_volumes,]
        
        if len(from_positions)==1:
            # ONE-TO-MANY
            from_positions = from_positions*len(to_positions)
        if len(to_positions)==1:
            # MANY-TO-ONE
            to_positions = to_positions*len(from_positions)
        if len(transfer_volumes)==1:
            # IF ONLY ONE VOLUME IS GIVEN, SET ALL VOLUMES EQUAL
            transfer_volumes = transfer_volumes*len(from_positions)
        if len(pre_transfer_delays)==1:
            # IF ONLY ONE DELAY IS GIVEN, SET ALL DELAY EQUAL
            pre_transfer_delays = pre_transfer_delays*len(from_positions)
            
            
        if repeat_each>1:
            from_positions = [x for x in from_positions for i in range(repeat_each)]
            transfer_volumes = [x for x in transfer_volumes for i in range(repeat_each)]

        from_positions = from_positions*repeat_all
        transfer_volumes = transfer_volumes*repeat_all
        
        num_transfers = min(len(from_positions),
                            len(to_positions),
                            len(transfer_volumes))
        
        return num_transfers,from_plate,from_positions,to_plate,to_positions,transfer_volumes,pre_transfer_delays
    
    
    def run(self, **kwargs):
        """ Run the transfer task.
    
        Kwargs:
            uncertainty_type (str). Should be 'gaussian' or 'uniform'.
            percentage_uncertainties (float)
            verbose (bool).  
            enforce_volume_limits (bool)
            volume_increment (float)   (optional. in nL)
            rotate_dest_180 (bool)   (optional, simulates a 180-degree rotation in the destination plate)
            
        Returns:
            None.
            
        Raises:
           Exception.
    
        Example use:
    
        XXXXXXXXXXXXXXXXXXXXXXX
    
        """         
        
        rotate_dest_180 = kwargs.get('rotate_dest_180',False)
        
        num_transfers,from_plate,from_positions,to_plate,to_positions,transfer_volumes,pre_transfer_delays = self.unpack_transfers(rotate_dest_180=rotate_dest_180)

        if len(to_positions)==0:
            print('\n** WARNING: running a TransferTask with zero transfers. **')
        
        if len(from_positions)>len(to_positions):
            print('\n** WARNING: running a TransferTask with more sources than destinations. Not all sources will be transferred. **')

        # round each transfer to the nearest increment as supported by the liquid handler
        volume_increment = kwargs.get('volume_increment',1)

        verbose = kwargs.get('verbose', True)
        uncertainties = kwargs.get('uncertainty_type', None)
        percentage_uncertainties = kwargs.get('percentage_uncertainties', None)
        enforce_volume_limits = kwargs.get('enforce_volume_limits',True)

        transfer_group_label = self._params.get('transfer_group_label',None)

        if uncertainties is not None and percentage_uncertainties is not None:
            transfer_volumes = np.array(transfer_volumes, dtype='float32')
            if uncertainties=='gaussian':
                transfer_volumes *= 1 + np.random.randn()*percentage_uncertainties
            elif uncertainties=='uniform':
                transfer_volumes *= 1 + np.random.uniform(-percentage_uncertainties,
                                                        percentage_uncertainties)
            else:
                raise Exception('Unknown distribution: %s' %(uncertainties,))
            transfer_volumes = list(np.int32(transfer_volumes))

        if(verbose):
            print('\nRunning transfer. len=%d' % (len(from_positions)) )
        
        for from_pos,to_pos,transfer_vol in zip(from_positions,to_positions,transfer_volumes):

            current_mixture_name = from_plate[from_pos]['mixture_name']
            source_min_volume = from_plate.get_param('volume_min')
            dest_max_volume = to_plate.get_param('volume_max')

            # enforce discrete volume increments
            round_vol = (transfer_vol // volume_increment) * volume_increment

            # check for available source volume
            source_avail_volume = from_plate[from_pos]['volume'] - source_min_volume
            if enforce_volume_limits and source_avail_volume < transfer_vol:
                raise Exception('Source Plate Exception:   %s, position %s %s, is depleted. Aborting protocol.'\
                                 % (from_plate._description,
                                    from_pos,
                                    Containers.position_to_lettergrid(from_pos)))

            # check for available destination volume
            dest_avail_volume = dest_max_volume - to_plate[to_pos]['volume']
            if enforce_volume_limits and source_avail_volume < transfer_vol:
                raise Exception('Destination Plate Exception:   %s, position %s %s, will overflow. Aborting protocol.'\
                                 % (to_plate._description,
                                    to_pos,
                                    Containers.position_to_lettergrid(to_pos)))
                
            # subtract volume from the source
            from_plate.add_compounds_to_location(
                        from_pos,
                        from_plate[from_pos]['compound_list'],
                        volume=-round_vol)
    
            # add volume to the destination
            to_plate.add_compounds_to_location(
                                to_pos,
                                from_plate[from_pos]['compound_list'],
                                volume=round_vol)
            if verbose:
                print('   transfer %.2E from %s to %s' % (round_vol,
                                                  str(from_pos),
                                                  str(to_pos))) 
            
            if transfer_group_label is not None:
                # write the new transfer group label
                to_plate[to_pos]['transfer group']=transfer_group_label
            else:
                pass
                # propagate the old label
                #to_plate[to_pos]['transfer group']=from_plate[from_pos].get('transfer group','')
    
    
    
        if(verbose):                        
            print('Transfer finished.\n')




        
# ============================================================

class TaskList:

    def __init__(self,*args,**kwargs):
        self._task_list = []
        self._description = kwargs.get('description','[no description]')
    
    def add(self,newtask):
        self._task_list.append(newtask)

    def merge(self,othertasklist):
        self._task_list.extend(othertasklist._task_list)

    def summarize(self):
        transfer_rate = 660.5; # Echo's empirical avgerage transfer rate [nL/s]
        total_transfers = 0
        total_volume = 0;
        from_plates = set()
        to_plates = set()
        for mytask in self._task_list:
            num_transfers,from_plate,from_positions,to_plate,to_positions,transfer_volumes,pre_transfer_delays = mytask.unpack_transfers()
            total_transfers += num_transfers
            total_volume += sum(transfer_volumes);
            from_plates.add(from_plate._description)
            to_plates.add(to_plate._description)
        print('\nTaskList summary:  %s' % (self._description))
        print('  %u subtasks, %u total transfers, \n  %u source plates %s, \n  %u destination plates %s, \n %u min. to complete, %u uL will be moved' \
              % (len(self._task_list),
                 total_transfers,
                 len(from_plates),
                 str(from_plates),
                 len(to_plates),
                 str(to_plates),
                 round(total_volume/transfer_rate/60),                                  
                 round(total_volume/1e3),                                                   
                 ))

    def run(self,**kwargs):
        #print('Run TaskList')
        for t in self._task_list:
            if t.has_param('description'):
                print('- Run Task:',t.get_param('description'))
            t.run(**kwargs)


# ------------------------------------------
class rowsum_TaskList(TaskList):
    
    def __init__(self,*args,**kwargs):
        TaskList.__init__(self,*args,**kwargs)
        
        from_plate=kwargs['from_plate']
        from_top_left=kwargs['from_top_left']
        from_bottom_right=kwargs['from_bottom_right']

        to_plate=kwargs['to_plate']
        to_top_left=kwargs['to_top_left']
        to_bottom_right=kwargs['to_bottom_right']
        
        transfer_volume=kwargs['transfer_volume']
        
        from_positions = Containers.bounds_to_positions( [ (from_top_left,from_bottom_right) ] )
        to_positions = Containers.bounds_to_positions( [ (to_top_left,to_bottom_right) ] )

        for i,thisrow in enumerate(Containers.separate_by_row(from_positions)):
            
            tt1 = TransferTask(from_plate=from_plate,
                                from_positions=thisrow,
                                to_plate=to_plate,
                                to_positions=to_positions[i],
                                transfer_volumes=transfer_volume,     # nL
                                transfer_group_label='rowsums'
                                )
            self.add(tt1)
            
# ------------------------------------------
class blocksum_TaskList(TaskList):
    
    def __init__(self,*args,**kwargs):
        TaskList.__init__(self,*args,**kwargs)
        
        from_plate=kwargs['from_plate']
        from_top_left=kwargs['from_top_left']
        from_bottom_right=kwargs['from_bottom_right']

        to_plate=kwargs['to_plate']
        to_position=kwargs['to_position']
        
        transfer_volume=kwargs['transfer_volume']
        
        from_positions = Containers.bounds_to_positions( [ (from_top_left,from_bottom_right) ] )

            
        tt1 = TransferTask(from_plate=from_plate,
                            from_positions=from_positions,
                            to_plate=to_plate,
                            to_positions=to_position,
                            transfer_volumes=transfer_volume,     # nL
                            transfer_group_label='blocksum'
                            )
        self.add(tt1)
            

# ------------------------------------------
class maldispot_TaskList(TaskList):
    
    def __init__(self,*args,**kwargs):
        TaskList.__init__(self,*args,**kwargs)
        
        from_plate=kwargs['from_plate']
        from_top_left=kwargs.get('from_top_left',None)          # OPTIONAL
        from_bottom_right=kwargs.get('from_bottom_right',None)  # OPTIONAL

        to_plate=kwargs['to_plate']
        to_top_left=kwargs.get('to_top_left',None)
        to_bottom_right=kwargs.get('to_bottom_right',None)
        
        transfer_volume=kwargs['transfer_volume']
        
        same_positions=kwargs.get('same_positions',False)
        repeat_each=kwargs.get('repeat_each',1)
        
        group_label=kwargs.get('group_label',None)
        
        
        if from_top_left is not None and from_bottom_right is not None:
            from_positions = Containers.bounds_to_positions( [ (from_top_left,from_bottom_right) ] )
        else:
            from_positions = from_plate.list_filled_positions()

        if same_positions is True:
            to_positions=from_positions
        else:
            if to_top_left is not None and to_bottom_right is not None:
                to_positions = Containers.bounds_to_positions( [ (to_top_left,to_bottom_right) ] )
            else:
                to_positions = to_plate.list_empty_positions()
            
        tt1 = TransferTask(from_plate=from_plate,
                            from_positions=from_positions,
                            to_plate=to_plate,
                            to_positions=to_positions,
                            transfer_volumes=transfer_volume,     # nL
                            transfer_group_label=group_label,
                            repeat_each=repeat_each
                            )
        self.add(tt1)
                        



# ------------------------------------------
class dilution_TaskList(TaskList):
    
    def __init__(self,*args,**kwargs):
        TaskList.__init__(self,*args,**kwargs)
        
        from_plate=kwargs['from_plate']
        from_positions=kwargs.get('from_positions')
        solvent_positions=kwargs.get('solvent_positions')

        to_plate=kwargs['to_plate']
        to_top_left=kwargs.get('to_top_left',None)          # optional
        to_bottom_right=kwargs.get('to_bottom_right',None)  # optional
        
        total_volume=kwargs['total_volume']
        dilution_ratios=kwargs['dilution_ratios']    # 0<x<1
                
        group_label=kwargs.get('group_label','dilution')
        

        # take from solvents multiple times if needed
        while len(solvent_positions)<(len(dilution_ratios)):
            solvent_positions = solvent_positions*2
        solvent_positions = solvent_positions[:len(dilution_ratios)]
        print('SOLVENTS',len(solvent_positions),solvent_positions)

        if to_top_left is not None and to_bottom_right is not None:
            to_positions = Containers.bounds_to_positions( [ (to_top_left,to_bottom_right) ] )
        else:
            to_positions = to_plate.list_empty_positions()

        # do a separate dilution series for each one
        for from_pos,to_pos in zip(from_positions,Containers.separate_by_count(to_positions,len(dilution_ratios))):
            tt1 = TransferTask(from_plate=from_plate,
                                from_positions=from_pos,
                                to_plate=to_plate,
                                to_positions=to_pos,
                                transfer_volumes=total_volume*dilution_ratios,     # nL
                                transfer_group_label= '%s (%s)' % (group_label,str(from_pos))
                                )
            self.add(tt1)
          
            tt1 = TransferTask(from_plate=from_plate,
                                from_positions=solvent_positions,
                                to_plate=to_plate,
                                to_positions=to_pos,
                                transfer_volumes=total_volume*(1-dilution_ratios),     # nL
                                transfer_group_label= '%s (%s)' % (group_label,str(from_pos)),
                                )
            self.add(tt1)
            
            
            
