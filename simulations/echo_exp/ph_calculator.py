import math

def calculate_resulting_ph_vol_conc(acid_vol, acid_conc, base_vol, base_conc): #L, molar
    acid_num_moles = acid_conc * acid_vol # moles
    base_num_moles = base_conc * base_vol # moles
    acid_is_limiting_agent = base_num_moles >= acid_num_moles
    if acid_is_limiting_agent:
        remaining_base_moles = base_num_moles - acid_num_moles
        oh_conc = remaining_base_moles / (acid_vol + base_vol)
        if remaining_base_moles == 0:
            h_conc = 1e-7
        else:
            h_conc = 1e-14 / oh_conc
    else:
        remaining_acid_moles = acid_num_moles - base_num_moles
        if remaining_base_moles == 0:
            h_conc = 1e-7
        else:
            h_conc = remaining_acid_moles / (acid_vol + base_vol)

    result_ph = -math.log10(h_conc)
    return result_ph

def calculate_acid_ph(acid_conc):
    return -math.log10(acid_conc)

def calculate_base_ph(base_conc):
    return 14.0 + math.log10(base_conc)

def calculate_similar_solutions_ph(first_sol_vol, first_sol_ph, second_sol_vol, second_sol_ph):
    first_is_water = abs(first_sol_ph - 7.0) < 0.9
    second_is_water = abs(second_sol_ph - 7.0) < 0.9
    if first_is_water and second_is_water:
        return 7.0
    elif first_is_water:
        return calculate_similar_solutions_ph(second_sol_vol, second_sol_ph, first_sol_vol, first_sol_ph)
    elif second_is_water:
        if first_sol_ph < 7:
            h_conc = 10**-first_sol_ph
            acid_num_moles = h_conc * first_sol_vol
            total_vol = first_sol_vol + second_sol_vol
            new_h_conc = acid_num_moles / total_vol
            return -math.log10(new_h_conc)
        else:
            h_conc = 10**-first_sol_ph
            oh_conc = 1e-14/h_conc
            base_num_moles = oh_conc * first_sol_vol
            total_vol = first_sol_vol + second_sol_vol
            new_oh_conc = base_num_moles / total_vol
            acid_num_moles = 1e-14/new_oh_conc
            ph = -math.log10(acid_num_moles)
            return ph
    else:
        first_h_conc = 10**-first_sol_ph
        second_h_conc = 10**-second_sol_ph
        first_num_moles = first_h_conc * first_sol_vol
        second_num_moles = second_h_conc * second_sol_vol
        total_num_moles = first_num_moles + second_num_moles
        h_combined_conc = total_num_moles / (first_sol_vol + second_sol_vol)
        ph = -math.log10(h_combined_conc)
        return ph





def calculate_resulting_ph_vol_ph(acid_vol, acid_ph, base_vol, base_ph): #L, ph
    if acid_ph > 7.0 and base_ph < 7.0:
        return calculate_resulting_ph_vol_ph(base_vol, base_ph, acid_vol, acid_ph)
    acid_h_conc = 10**-acid_ph
    acid_num_moles = acid_h_conc * acid_vol
    base_h_conc = 10**-base_ph
    oh_conc = 1e-14/(base_h_conc)
    base_num_moles = oh_conc * base_vol
    acid_is_limiting_agent = base_num_moles > acid_num_moles
    if (abs(base_num_moles - acid_num_moles)/max(base_num_moles, acid_num_moles) < 1e-8):
        return 7.00
    if acid_is_limiting_agent:
        remaining_base_moles = base_num_moles - acid_num_moles
        oh_conc = remaining_base_moles / (acid_vol + base_vol)
        h_conc = 1e-14 / (oh_conc)
    else:
        remaining_acid_moles = acid_num_moles - base_num_moles
        h_conc = remaining_acid_moles / (acid_vol + base_vol)
    result_ph = -math.log10(h_conc)
    return result_ph


def calculate_buffer_ph(ka, acid_conc, base_conc):
    pka = -math.log10(ka)
    return (pka + math.log10(base_conc / acid_conc)), acid_conc, base_conc


def calculate_buffer_with_base(ka, buffer_acid_conc, buffer_base_conc, base_conc):
    new_acid_conc = buffer_acid_conc - base_conc
    new_base_conc = buffer_base_conc + base_conc
    return calculate_buffer_ph(ka, new_acid_conc, new_base_conc)

def calculate_buffer_with_acid(ka, buffer_acid_conc, buffer_base_conc, acid_conc):
    new_acid_conc = buffer_acid_conc + acid_conc
    new_base_conc = buffer_base_conc - acid_conc
    return calculate_buffer_ph(ka, new_acid_conc, new_base_conc)

def calculate_buffer_with_buffer(ka, first_buffer_acid_conc, first_buffer_base_conc, second_buffer_acid_conc, second_buffer_base_conc):
    new_acid_conc = first_buffer_acid_conc + second_buffer_acid_conc
    new_base_conc = first_buffer_base_conc + second_buffer_base_conc
    return calculate_buffer_ph(ka, new_acid_conc, new_base_conc)