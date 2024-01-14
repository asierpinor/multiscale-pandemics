"""
This module contains a collection of functions to load USA data on population, commuting and airtravel,
apply some data manipulations to make the data mutually compatible and use this data to compute a matrix T
of transmission probabilities.
"""

import numpy as np
from scipy.sparse import csr_array
import pandas as pd

# Default parameters
DEFAULT_BETA = 1
DEFAULT_AIR_PREFACTOR = 1
DEFAULT_COMM_PREFACTOR = 1/4

# Data folders
data_folder_FAA = '../data/FAA_data'
data_folder_census = '../data/UScensus_data'
data_folder_simplemaps = '../data/Simplemaps_data'
data_folder_statesFIPS = '../data/StatesFIPS_data'

# Name of data files
census_year = 2021
census_file = 'population-bycounty-%i-UScensus.csv'%(census_year)

commuting_year = '2011-15'
commuting_file = 'commuting_countyLevel_sortResidence_UScensus_edited.xlsx'

USstates_file = 'USstates_FIPS.xlsx'

UScities_file = 'uscities_to_county.xlsx'

faa_year = '19'
faa_airtravel_file = 'cy%s-commercial-service-enplanements_edited.xlsx'%(faa_year)



def load_and_process_usa_census_data (data_folder=data_folder_census, data_file=census_file):
    """Load census data and separate into counties and states."""

    filename = '%s/%s'%(data_folder,data_file)

    # Read data
    census = pd.read_csv(filename,engine='python',encoding='latin1',comment='#',dtype='string')

    # Keep first columns only
    census = census.drop(columns=census.columns[list(range(10,35))])
    census = census.drop(columns=census.columns[[7,8]]) # Keep only 2021 estimate, column 9

    # Change some column names
    census = census.rename(columns={'POPESTIMATE2021': 'POPULATION'})

    # Extract just counties/states
    states_census = census[census.iloc[:,0] == '040']
    counties_census = census[census.iloc[:,0] == '050']

    # Remove county/states column
    states_census = states_census.drop(columns=['SUMLEV'])
    counties_census = counties_census.drop(columns=['SUMLEV'])

    # Add FIPS column
    counties_census['FIPS'] = counties_census['STATE'] + counties_census['COUNTY']
    states_census['FIPS'] = states_census['STATE']

    return census, counties_census, states_census



def load_and_process_commuting_data (data_folder=data_folder_census, data_file=commuting_file):
    """Load commuting data, remove some entries, and simplify."""

    filename = '%s/%s'%(data_folder, data_file)

    # Read data
    commuting = pd.read_excel(filename, engine='openpyxl', dtype=str)

    # Remove NaNs related to rows indicating travel outside of the US
    commuting = commuting[commuting.iloc[:,4].notnull()]
    commuting = commuting[commuting.iloc[:,5].notnull()]

    # Remove Puerto Rico data
    commuting = commuting[commuting.iloc[:,2] != 'Puerto Rico']
    commuting = commuting[commuting.iloc[:,6] != 'Puerto Rico']

    # Remove redundant columns
    commuting = commuting.drop(commuting.columns[[2,3,6,7,9]], axis=1)

    # Rename columns
    commuting.columns = ['State FIPS res', 'County FIPS res', 'State FIPS work', 'County FIPS work', 'Flow']

    # Combine FIPS codes into one column
    temp_FIPS_res = commuting.iloc[:,0]+commuting.iloc[:,1]
    temp_FIPS_work = (commuting.iloc[:,2]+commuting.iloc[:,3])  # Replaces occurrences of "0" n=1 times
    temp_FIPS_work = temp_FIPS_work.str.replace('0', '',n=1)

    commuting["FIPS residency"] = temp_FIPS_res
    commuting["FIPS work"] = temp_FIPS_work

    # Keep only 3 columns
    commuting = commuting[['FIPS residency','FIPS work','Flow']]

    # Change type
    commuting['Flow'] = commuting['Flow'].astype(float)

    return commuting



def rename_old_counties_commuting (commuting, counties_census):
    """Update Valdez commuting data to current Chugach and Copper."""

    # Copy Valdez commuting data
    comm_Valdez = pd.concat([ commuting[commuting['FIPS work'] == '02261'] ,
                              commuting[commuting['FIPS residency'] == '02261'] ], ignore_index=True)

    comm_Chugach = comm_Valdez.copy()
    comm_Copper = comm_Valdez.copy()

    # Extract population of Chugach and Copper
    pop_Chugach = int( counties_census[ counties_census['FIPS'] == '02063' ]['POPULATION'].iloc[0] )
    pop_Copper = int( counties_census[ counties_census['FIPS'] == '02066' ]['POPULATION'].iloc[0] )
    pop_CC = pop_Chugach + pop_Copper

    # Chugach coommuting data from Valdez with flow rescaled by population
    Valdez_to_Chugach = {'02261': '02063'}
    comm_Chugach['Flow'] = comm_Chugach['Flow']*pop_Chugach/pop_CC
    comm_Chugach['FIPS residency'] = comm_Chugach['FIPS residency'].replace(Valdez_to_Chugach)  # inplace=True
    comm_Chugach['FIPS work'] = comm_Chugach['FIPS work'].replace(Valdez_to_Chugach)   # inplace=True

    # Copper coommuting data from Valdez with flow rescaled by population
    Valdez_to_Copper = {'02261': '02066'}
    comm_Copper['Flow'] = comm_Copper['Flow']*pop_Copper/pop_CC
    comm_Copper['FIPS residency'] = comm_Copper['FIPS residency'].replace(Valdez_to_Copper)  # inplace=True
    comm_Copper['FIPS work'] = comm_Copper['FIPS work'].replace(Valdez_to_Copper)   # inplace=True

    ### NOTE: There will be no commuting between Copper and Chugach. I could compute something from the value of
    ###       self-commuting, but this is completely negligible

    # Combine Chuggar and Copper
    comm_ChugachCopper = pd.concat( [ comm_Chugach , comm_Copper ] , ignore_index=True)

    # Remove Valdez from commuting file
    commuting = commuting[ commuting['FIPS work'] != '02261' ]
    commuting = commuting[ commuting['FIPS residency'] != '02261' ]

    # Add Chuggar and Copper
    commuting = pd.concat( [commuting, comm_ChugachCopper] , ignore_index=True )

    return commuting



def remove_nonmatching_FIPS_census_commuting (commuting, counties_census):
    """Remove counties that don't appear in both census and commuting (if any)."""
    
    # Extract counties FIPS codes from census file
    pop_countiesFIPS = (counties_census.loc[:,'STATE']+counties_census.loc[:,'COUNTY']).to_numpy(dtype=str)
    
    # Create list of FIPS keys present in this file to crosscheck with FIPS keys of population file
    commuting_countiesFIPS = []

    for ii in range(len(commuting["FIPS residency"])):
        if commuting["FIPS residency"].iloc[ii] not in commuting_countiesFIPS:
            commuting_countiesFIPS.append(commuting["FIPS residency"].iloc[ii])
        if commuting["FIPS work"].iloc[ii] not in commuting_countiesFIPS:
            commuting_countiesFIPS.append(commuting["FIPS work"].iloc[ii])
    
    #### Check inconsistencies between two lists ####

    remove_FIPS_from_population = []
    #print("FIPS in population file NOT in commuting file")
    for key in pop_countiesFIPS:
        if key not in commuting_countiesFIPS:
            remove_FIPS_from_population.append(key)
    #print(remove_FIPS_from_population)

    remove_FIPS_from_commuting = []
    #print("FIPS in commuting file NOT in population file")
    for key in commuting_countiesFIPS:
        if key not in pop_countiesFIPS:
            remove_FIPS_from_commuting.append(key)
    #print(remove_FIPS_from_commuting)

    #### Remove entries from both arrays ####

    select_pop_rows = [ pop_countiesFIPS[ii] not in remove_FIPS_from_population for ii in range(len(pop_countiesFIPS)) ]
    counties_census = counties_census.iloc[ select_pop_rows , :]

    commuting = commuting.iloc[ [ commuting.iloc[ii,0] not in remove_FIPS_from_commuting for ii in range(len(commuting)) ] , :]
    commuting = commuting.iloc[ [ commuting.iloc[ii,1] not in remove_FIPS_from_commuting for ii in range(len(commuting)) ] , :]
    
    return commuting, counties_census
    



def set_FIPS (counties_census, states_census):
    """Define FIPS arrays and index-FIPS dictionaries for counties, states and divisions."""

    #### Define final FIPS arrays ####

    counties_FIPS = (counties_census.loc[:,'STATE']+counties_census.loc[:,'COUNTY']).to_numpy(dtype=str)
    states_FIPS = (states_census.loc[:,'STATE']).to_numpy(dtype=str)
    divisions_FIPS = np.array( np.arange(1,10) , dtype=str)
            
    #### Define dictionaries FIPS to index ####

    countyFIPS_to_index = { counties_FIPS[ii]: ii for ii in range(len(counties_FIPS)) }
    statesFIPS_to_index = { states_FIPS[ii]: ii for ii in range(len(states_FIPS)) }
    divisionsFIPS_to_index = { divisions_FIPS[ii]: ii for ii in range(len(divisions_FIPS)) }

    index_to_countyFIPS = { ii: counties_FIPS[ii] for ii in range(len(counties_FIPS)) }
    index_to_statesFIPS = { ii: states_FIPS[ii] for ii in range(len(states_FIPS)) }
    index_to_divisionsFIPS = { ii: divisions_FIPS[ii] for ii in range(len(divisions_FIPS)) }

    return (counties_FIPS, states_FIPS, divisions_FIPS, countyFIPS_to_index, statesFIPS_to_index, divisionsFIPS_to_index,
            index_to_countyFIPS, index_to_statesFIPS, index_to_divisionsFIPS)


def load_names_states_divisions (divisions_FIPS, data_folder=data_folder_statesFIPS, data_file=USstates_file):
    """Load US states acronyms and define US divisions names."""

    filename = '%s/%s'%(data_folder, data_file)

    USstates = pd.read_excel(filename, engine='openpyxl', dtype=str)

    # Remove Puerto Rico and Virgin Islands
    USstates = USstates[ USstates['ST-Code'] != 'PR' ]
    USstates = USstates[ USstates['ST-Code'] != 'VI' ]

    USdivisions = pd.DataFrame()
    USdivisions['FIPS'] = divisions_FIPS
    div_names = { '1': 'New England', '2': 'Middle Atlantic', '3': 'East North Central', '4': 'West North Central', '5': 'South Atlantic',\
                '6': 'East South Central', '7': 'West South Central', '8': 'Mountain', '9': 'Pacific' }
    USdivisions['Division'] = USdivisions['FIPS'].replace(div_names)

    return USstates, USdivisions


def define_population_arrays (counties_census, states_census, USdivisions):
    """
    Define population arrays of counties, states and divisions from data,
    and define indices relating regions to sub-regions.
    """

    # Define just population arrays
    states_pop = states_census.loc[:,'POPULATION'].to_numpy(dtype=int)
    counties_pop = counties_census.loc[:,'POPULATION'].to_numpy(dtype=int)

    # Combine state populations into divisions
    divisions_pop = []
    for nn in range(1,10):
        states_of_division_nn = states_census[states_census.loc[:,'DIVISION'] == str(nn)]
        states_of_division_nn = states_of_division_nn.loc[:,'POPULATION'].to_numpy(dtype=int)
        
        divisions_pop.append(np.sum(states_of_division_nn))
        
    divisions_pop = np.array(divisions_pop)

    # Indices of counties within state
    counties_in_state = {}
    for ii in range(len(states_pop)):
        temp_counties_inds = (counties_census.loc[:,'STATE'] == states_census.loc[states_census.index.values[ii],'STATE']).to_numpy(dtype=bool)
        counties_in_state[ii] = np.arange(len(counties_pop))[ temp_counties_inds ]
        
    # Indices of states within divisions
    states_in_division = {}
    for ii in range(len(divisions_pop)):
        temp_states_inds = (states_census.loc[:,'DIVISION'] == USdivisions.loc[USdivisions.index.values[ii],'FIPS']).to_numpy(dtype=bool)
        states_in_division[ii] = np.arange(len(states_pop))[ temp_states_inds ]

    return counties_pop, states_pop, divisions_pop, counties_in_state, states_in_division



def transform_commuting_into_matrix (commuting, n_counties, countyFIPS_to_index):
    """Transform commuting data into matrix form."""

    # Transform FIPS into county index
    commuting['Index FIPS residency'] = commuting['FIPS residency'].replace(countyFIPS_to_index)  # inplace=True
    commuting['Index FIPS work'] = commuting['FIPS work'].replace(countyFIPS_to_index)   # inplace=True

    # Transform commuting into sparse matrix
    commuting_matrix = csr_array((commuting['Flow'].to_numpy(dtype=float), \
                (commuting['Index FIPS work'].to_numpy(dtype=int), commuting['Index FIPS residency'].to_numpy(dtype=int))),\
            shape=(n_counties, n_counties))
    # Format: data, rows, columns

    # Remove diagonal (self-commuting)
    for ii in range(commuting_matrix.shape[0]):
        commuting_matrix[ii,ii] = 0

    return commuting_matrix


def load_UScity_to_county (data_folder=data_folder_simplemaps, data_file=UScities_file):
    """Load US cities data and define city-to-county dictionary."""

    filename = '%s/%s'%(data_folder, data_file)

    uscities = pd.read_excel(filename, engine='openpyxl', dtype=str)

    uscities = uscities.drop(uscities.columns[range(6,17)], axis=1) # Use if column labels are not numbers
    uscities = uscities.drop(uscities.columns[1], axis=1)

    # Add '0' to the beginning of FIPS that are only 4 chars long (due to numerical value being < 10000)
    missing_zero_fips = [len(uscities['county_fips'].iloc[ii])==4 for ii in range(len(uscities))]
    uscities['county_fips'].iloc[ missing_zero_fips ] = '0' + uscities['county_fips'].iloc[ missing_zero_fips ]

    citytocounty = { uscities['city'].iloc[ii]+'-'+uscities['state_id'].iloc[ii]: uscities['county_fips'].iloc[ii] for ii in range(len(uscities)) }

    return citytocounty


def load_and_process_flight_data (USstates, citytocounty, counties_FIPS, countyFIPS_to_index, n_counties,
                                  data_folder=data_folder_FAA, data_file=faa_airtravel_file):
    """Load and edit enplanements data for each airport and create array of enplanements per county."""

    filename = '%s/%s'%(data_folder, data_file)

    faatravel = pd.read_excel(filename, engine='openpyxl', dtype=str)

    # Remove some columns and null rows
    faatravel = faatravel.drop(faatravel.columns[[0,1,9,10]], axis=1) # Use if column labels are not numbers
    faatravel = faatravel.drop([526,527])
    faatravel = faatravel[faatravel.loc[:,'CY 19 Enplanements'].notnull()]

    # Remove airports not in the 51 US states
    rows_in_US = [ faatravel['ST'].iloc[ii] in USstates['ST-Code'].to_numpy(dtype=str) for ii in range(len(faatravel)) ]
    faatravel = faatravel.iloc[ rows_in_US , : ]

    # Replace some small cities by big cities
    toreplace = ['Dulles'    ,'Windsor Locks','Kailua Kona','Grand Canyon','Bristol/Johnson/Kingsport','Bloomington-Normal Airport','Clinton (Township of)','Nantucket','Deadhorse'  ,'Barrow'   ,'Westerly' ,'Block Island',"St. Mary's (ANV/ANVSA)",'Provincetown','Bar Harbor','Eastsound','Devils Lake','Fort Leonard Wood (U.S. Army)','Wade Hampton (Census Area)','Barrow (County)','Kalskag'      ,'St. Michael (ANV/ANVSA)','Adak (Naval) Station/Mitchell Field','Saint Paul Island','Islip'        ,'Newburgh' ,'Hyannis'      ,'Massena']
    replaceby = ['Washington','Hartford'     ,'Waimea'     ,'Williams'    ,'Kingsport'                ,'Bloomington'               ,'Lansing'              ,'Madaket'  ,'Prudhoe Bay','Utqiagvik','Weekapaug','Weekapaug'   ,'Russian Mission'       ,'Harwich Port','Ellsworth' ,'Anacortes','Fort Totten','Fort Leonard Wood'            ,'Russian Mission'           ,'Wainwright'     ,'Upper Kalskag','St. Michael'            ,'Adak'                               ,'St. Paul'         ,'Central Islip','Balmville','West Yarmouth','Brasher Falls']
    faatravel['City'] = faatravel['City'].replace(toreplace,replaceby)

    # Replace Washington-VA by Washington DC
    faatravel.loc[ faatravel['City'] == 'Washington', 'ST' ] = 'DC'

    # Change to city-state, so cities can be uniquely identified in citytocounty dict.
    faatravel['City'] = faatravel['City'] + '-' + faatravel['ST']

    # Add FIPS column by mapping city to FIPS
    faatravel['FIPS'] = faatravel['City'].replace(citytocounty)

    # Rename column
    faatravel = faatravel.rename(columns={'CY 19 Enplanements': 'Enplanements'})

    # Change data type to int
    faatravel = faatravel.astype({'Enplanements': 'int'})

    # Remove airports whose FIPS does not coincide with list of county FIPS
    faatravel = faatravel[ faatravel['FIPS'].isin(counties_FIPS) ]

    # Make copy of faatravel
    #faatravel_full = faatravel

    # Keep only 2 columns
    faatravel = faatravel[['FIPS','Enplanements']]

    # Check if any FIPS is repeated and combine enplanements within same county
    row = 0
    while row<len(faatravel):
        
        ii = faatravel.index[row]
        
        check_fips = faatravel['FIPS'].loc[ii]
        ind_same_fips = (faatravel.index[faatravel['FIPS'] == check_fips]).to_numpy(dtype=int)
        ind_same_fips = ind_same_fips[ind_same_fips != ii] # Remove self
        
        add_enplanements = np.sum( [ faatravel['Enplanements'].loc[rr] for rr in ind_same_fips ] )
        faatravel.loc[ii,'Enplanements'] = faatravel.loc[ii,'Enplanements'] + int(add_enplanements)
        
        faatravel = faatravel.drop(ind_same_fips)
        
        row = row+1

    # Yearly to daily enplanements
    faatravel['Enplanements'] = (faatravel['Enplanements']/365) #.astype(int)   # Takes floor function by int conversion

    # Change FIPS to index
    faatravel['FIPS'].replace(countyFIPS_to_index, inplace=True)
    faatravel = faatravel.astype({'FIPS': 'int'})
    faatravel = faatravel.rename(columns={'FIPS': 'County'})

    # Transform into array
    airport_enplane = np.zeros(n_counties, dtype=int)
    for ii in range(len(faatravel)):
        airport_enplane[ faatravel['County'].iloc[ii] ] = faatravel['Enplanements'].iloc[ii]

    return faatravel, airport_enplane


def transform_airtravel_to_matrix (n_counties, airport_enplane):
    """Transform enplanements into flight matrix."""

    Ftravelers_matrix = np.zeros((n_counties,n_counties))

    total_enplane = np.sum(airport_enplane) # Number of enplanements per day

    # Define Fij matrix as number of people traveling from j to i. STEM method.
    for ii in range(n_counties):
        Ftravelers_matrix[ii,:] = airport_enplane * airport_enplane[ii] / total_enplane
        Ftravelers_matrix[ii,ii] = 0

    # Actual number of enplanements according to Fmatrix, and how many airports have incoming passengers
    mod_airport_enplane = np.squeeze(np.sum(Ftravelers_matrix,axis=0))
    mod_airport_arrivals = np.squeeze(np.sum(Ftravelers_matrix,axis=1))

    return Ftravelers_matrix, mod_airport_enplane, mod_airport_arrivals


def compute_mobility_matrix (n_counties, Ftravelers_matrix, commuting_matrix, counties_pop,
                             air_prefactor=DEFAULT_AIR_PREFACTOR, comm_prefactor=DEFAULT_COMM_PREFACTOR):
    """Combine commuting and airtravel into mobility matrix."""

    mobility_matrix = np.zeros((n_counties,n_counties))

    #t_disease = 4  # 7?
    #t_travel = 7  # 10?
    #air_prefactor = 1/4 * min(t_travel,t_disease)      # 1/2 due to transfers, 1/2 due to outbound
    #comm_prefactor = 1/2      # fraction of time a commuter spends outside county

    mobility_matrix = air_prefactor * Ftravelers_matrix + comm_prefactor * commuting_matrix

    # Fraction of mobile people
    frac_mobilers = np.squeeze(np.sum(mobility_matrix,axis=0)) / counties_pop

    # Compute actual number of people in each county at any moment in time
    counties_realpop = counties_pop - np.sum(mobility_matrix,axis=0) + np.sum(mobility_matrix,axis=1)
    frac_pop = counties_realpop/counties_pop

    # Set diagonal such that sum_i M_ij = N_j
    mobility_matrix += np.diag(counties_pop - np.sum(mobility_matrix,axis=0))

    return mobility_matrix, counties_realpop, frac_mobilers, frac_pop


def compute_Tij_from_mobility (n_counties, counties_realpop, counties_pop, mobility_matrix, beta=DEFAULT_BETA):
    """Compute transmission probability matrix T from mobility matrix."""

    T0_total = np.zeros((n_counties,n_counties))

    Plocal = beta / counties_realpop     # IMPORTANT: "realpop", not "pop"
    Piinj = mobility_matrix / counties_pop.reshape(1,len(counties_pop))  # Probability of i to be in county c_j

    # T0
    T0_total = np.transpose(Piinj) @ (Plocal.reshape(len(Plocal),1) * Piinj )

    return T0_total


def compute_Tij_from_data (output='basic',beta=DEFAULT_BETA, air_prefactor=DEFAULT_AIR_PREFACTOR, comm_prefactor=DEFAULT_COMM_PREFACTOR):
    """
    Load data, manipulate data, and compute transmission probability matrix T.
    
    This is a master function that sequentially calls all of the above functions.
    """

    census, counties_census, states_census = load_and_process_usa_census_data()

    commuting = load_and_process_commuting_data()
    commuting = rename_old_counties_commuting (commuting, counties_census)
    commuting, counties_census = remove_nonmatching_FIPS_census_commuting (commuting, counties_census)

    (counties_FIPS, states_FIPS, divisions_FIPS, countyFIPS_to_index,
    statesFIPS_to_index, divisionsFIPS_to_index,index_to_countyFIPS,
    index_to_statesFIPS, index_to_divisionsFIPS) = set_FIPS (counties_census, states_census)

    USstates, USdivisions = load_names_states_divisions (divisions_FIPS)

    (counties_pop, states_pop, divisions_pop,
    counties_in_state, states_in_division) = define_population_arrays (counties_census, states_census, USdivisions)
    n_counties = len(counties_pop)

    commuting_matrix = transform_commuting_into_matrix (commuting, n_counties, countyFIPS_to_index)

    citytocounty = load_UScity_to_county()

    faatravel, airport_enplane = load_and_process_flight_data (USstates, citytocounty, counties_FIPS,
                                                                    countyFIPS_to_index, n_counties)
    
    Ftravelers_matrix, mod_airport_enplane, mod_airport_arrivals = transform_airtravel_to_matrix (n_counties, airport_enplane)

    (mobility_matrix, counties_realpop,
    frac_mobilers, frac_pop) = compute_mobility_matrix (n_counties, Ftravelers_matrix, commuting_matrix, counties_pop,
                                                        air_prefactor=air_prefactor, comm_prefactor=comm_prefactor)
    T0_total = compute_Tij_from_mobility (n_counties, counties_realpop, counties_pop, mobility_matrix, beta=beta)

    if output == 'T0':
        return T0_total
    
    if output == 'basic':
        return T0_total, counties_pop, counties_in_state, states_in_division
    
    if output == 'all':
        return (census, counties_census, states_census, commuting, counties_FIPS, states_FIPS, divisions_FIPS, countyFIPS_to_index,
                statesFIPS_to_index, divisionsFIPS_to_index, index_to_countyFIPS, index_to_statesFIPS, index_to_divisionsFIPS,
                USstates, USdivisions, counties_pop, states_pop, divisions_pop, counties_in_state, states_in_division,
                n_counties, commuting_matrix, citytocounty, faatravel, airport_enplane, Ftravelers_matrix, mod_airport_enplane,
                mod_airport_arrivals, mobility_matrix, counties_realpop, frac_mobilers, frac_pop, T0_total)







