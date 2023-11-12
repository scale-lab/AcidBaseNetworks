import numpy as np

from ph_calculator import *

def read_file(filepath):
    with open(filepath, "r") as f:
        return f.read()

class AcidBaseNetwork:
    """" Describes the structure of an acid base network """ 
    def __init__(self, num_neurons=2, weights_per_neuron=64, img_width=8, img_height=8):
        self.num_neurons = num_neurons
        self.weights_per_neuron = weights_per_neuron
        self.neurons = [[] for _ in range(0, self.num_neurons)]
        self.img_width = img_width
        self.img_height = img_height
        self.img = None
        self.weights_loaded = False
        self.image_loaded = False
        self.neurons_outputs_acid_vol = []
        self.neurons_outputs_base_vol = []
        self.neurons_outputs_water_vol = []
        self.grayscale = False
        for i in range(0, self.num_neurons):
            self.neurons_outputs_acid_vol.append([0, 0])
            self.neurons_outputs_base_vol.append([0, 0])
            self.neurons_outputs_water_vol.append([0, 0])

    def load_weights(self, filename):
        """ Loads weights into neurons """
        weights = [int(float(val)) for val in read_file(filename).splitlines()]
        for i in range(0, len(self.neurons)):
            self.neurons[i] = weights[i * self.weights_per_neuron: (i + 1) * self.weights_per_neuron]
        self.weights_loaded = True

    
    def load_image(self, filename, img_levels=None):
        """ load an image file into the network """
        self.img = [int(float(val)) for val in read_file(filename).splitlines()]
        if img_levels is not None:
            # self.img = list(map(lambda el: img_levels.index(el) - 1,self.img))
            self.grayscale = True
        self.image_loaded = True

    def generate_data_plate(self, acid_source, base_source, starting_letter, starting_index, max_row, max_col):
        pos_list_data_source = [];
        pos_list_data_destin = [];
        letter_counter = ord(starting_letter)
        index_counter = starting_index
        get_next_counter = lambda letter, index: (letter, index + 1) if index + 1 <= max_col else (letter + 1, 1)
        for _ in range(0, self.num_neurons):
            for val in self.img:
                is_acid = val == 1
                first_dest_pos = chr(letter_counter) + str(index_counter)
                letter_counter, index_counter = get_next_counter(letter_counter, index_counter)
                second_dest_pos = chr(letter_counter) + str(index_counter)
                letter_counter, index_counter = get_next_counter(letter_counter, index_counter)
                if is_acid:
                    pos_list_data_source.append(acid_source)
                    pos_list_data_source.append(base_source)
                else:
                    pos_list_data_source.append(base_source)
                    pos_list_data_source.append(acid_source)
                pos_list_data_destin.append(first_dest_pos)
                pos_list_data_destin.append(second_dest_pos)
            letter_counter = letter_counter + 2
            index_counter = 1
        return pos_list_data_source, pos_list_data_destin

    def generate_weighted_data_plate(self, acid_positions, base_positions, water_positions, source_max_transfers, starting_letter, starting_index, max_row, max_col, transfer_unit_vol, acid_conc, base_conc):
        pos_list_data_source = [];
        pos_list_data_destin = [];
        letter_counter = ord(starting_letter)
        index_counter = starting_index
        get_next_counter = lambda letter, index: (letter, index + 1) if index + 1 <= max_col else (letter + 1, 1)
        weight_index = 0

        current_acid_index = 0
        current_base_index = 0
        current_water_index = 0
        base_ph = calculate_base_ph(base_conc)
        acid_ph = calculate_acid_ph(acid_conc)
        water_pos_list_data_source = []
        water_pos_list_data_destin = []
        for neuron_index in range(0, self.num_neurons):
            for val in self.img:
                flip = self.neurons[neuron_index][weight_index] == -1
                # flip = False
                weight_index = (weight_index + 1) % self.weights_per_neuron
                is_on = val > 0
                first_dest_pos = chr(letter_counter) + str(index_counter)
                letter_counter, index_counter = get_next_counter(letter_counter, index_counter)
                second_dest_pos = chr(letter_counter) + str(index_counter)
                letter_counter, index_counter = get_next_counter(letter_counter, index_counter)

                added_vol = transfer_unit_vol
                
                if self.grayscale and is_on:
                    added_vol = transfer_unit_vol * (abs(val)/4.0)


                if (is_on and not flip) or (not is_on and flip):
                    pos_list_data_source.append(acid_positions[current_acid_index])
                    pos_list_data_source.append(base_positions[current_base_index])
                    self.neurons_outputs_acid_vol[neuron_index][0] = self.neurons_outputs_acid_vol[neuron_index][0] + added_vol #left
                    self.neurons_outputs_base_vol[neuron_index][1] = self.neurons_outputs_base_vol[neuron_index][1] + added_vol #right
                else:
                    pos_list_data_source.append(base_positions[current_base_index])
                    pos_list_data_source.append(acid_positions[current_acid_index])
                    self.neurons_outputs_base_vol[neuron_index][0] = self.neurons_outputs_base_vol[neuron_index][0] + added_vol #left
                    self.neurons_outputs_acid_vol[neuron_index][1] = self.neurons_outputs_acid_vol[neuron_index][1] + added_vol #right
                


                pos_list_data_destin.append(first_dest_pos)
                pos_list_data_destin.append(second_dest_pos)

                if self.grayscale and is_on:
                    for _ in range(0, val):
                        water_pos_list_data_source.append(water_positions[current_water_index])
                        water_pos_list_data_destin.append(first_dest_pos)
                        
                        current_water_index = (current_water_index + 1) % len(water_positions)
                        
                        water_pos_list_data_source.append(water_positions[current_water_index])
                        water_pos_list_data_destin.append(second_dest_pos)

                        current_water_index = (current_water_index + 1) % len(water_positions)

                        self.neurons_outputs_water_vol[neuron_index][0] = self.neurons_outputs_water_vol[neuron_index][0] + transfer_unit_vol
                        self.neurons_outputs_water_vol[neuron_index][1] = self.neurons_outputs_water_vol[neuron_index][1] + transfer_unit_vol


                current_acid_index = (current_acid_index + 1) % len(acid_positions)
                current_base_index = (current_base_index + 1) % len(base_positions)



            letter_counter = letter_counter + 2
            index_counter = 1
        return pos_list_data_source, pos_list_data_destin
            
    def generate_weights_plate(self, source_image, starting_letter, starting_index, max_row, max_col):
        pos_list_data_source = [];
        pos_list_data_destin = [];
        letter_counter = ord(starting_letter)
        index_counter = starting_index
        get_next_counter = lambda letter, index: (letter, index + 1) if index + 1 <= max_col else (letter + 1, 1)
        i = 0
        weight_index = 0
        while i < len(source_image):
            neuron_index = int(i / (2 * self.weights_per_neuron))
            flip = self.neurons[neuron_index][weight_index] == -1
            weight_index = (weight_index + 1) % self.weights_per_neuron

            first_dest_pos = chr(letter_counter) + str(index_counter)
            letter_counter, index_counter = get_next_counter(letter_counter, index_counter)
            second_dest_pos = chr(letter_counter) + str(index_counter)
            letter_counter, index_counter = get_next_counter(letter_counter, index_counter)
            if flip:
                pos_list_data_source.append(source_image[i + 1])
                pos_list_data_source.append(source_image[i])
            else:
                pos_list_data_source.append(source_image[i])
                pos_list_data_source.append(source_image[i + 1])
            pos_list_data_destin.append(first_dest_pos)
            pos_list_data_destin.append(second_dest_pos)
            
            if i > 0 and ((i + 2) % (self.img_width * self.img_height * 2) == 0):
                letter_counter = letter_counter + 2
                index_counter = 1
            i = i + 2
        return pos_list_data_source, pos_list_data_destin

    def generate_summation_plate(self, source_image, starting_letter, starting_index, max_row, max_col):
        pos_list_data_source = [];
        pos_list_data_destin = [];
        mixed = {}
        get_next_counter = lambda letter, index: (letter, index + 1) if index + 1 <= max_col else (letter + 1, 1)
        i = 0
        while i < len(source_image):
            neuron_index = int(i / (2 * self.weights_per_neuron))
            first_dest_pos = chr(ord(starting_letter) + neuron_index) + str(starting_index)
            second_dest_pos = chr(ord(starting_letter) + neuron_index) + str(starting_index + 1)
            pos_list_data_source.append(source_image[i])
            pos_list_data_source.append(source_image[i + 1])
            pos_list_data_destin.append(first_dest_pos)
            pos_list_data_destin.append(second_dest_pos)
            
            i = i + 2
        return pos_list_data_source, pos_list_data_destin

    def get_expected_outputs(self, acid_ph=1, base_ph=13):
        result = []
        for i in range(0, self.num_neurons):
            left = round(calculate_resulting_ph_vol_ph(self.neurons_outputs_acid_vol[i][0], acid_ph, self.neurons_outputs_base_vol[i][0], base_ph), 2)
            right = round(calculate_resulting_ph_vol_ph(self.neurons_outputs_acid_vol[i][1], acid_ph, self.neurons_outputs_base_vol[i][1], base_ph), 2)
            result.append([left, right])
        return result

    def generate_pool_plate(self, source_image, indicator_source, max_row, max_col):
        pos_list_data_source = [];
        pos_list_data_destin = [];
        for well in source_image:
            pos_list_data_source.append(indicator_source)
            pos_list_data_destin.append(well)
        return pos_list_data_source, pos_list_data_destin 


            
            


