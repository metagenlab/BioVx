#! /usr/bin/env python

import pickle


'''
The purpose of this program is to update the dictionary which contains the information regarding substrates of several proteins.
'''


def update_dict(fileNumber):


	#Open the file 'tcdb_substrates.tsv' and advance beyond the first line, which contains the header information
    
        fileName ='ChEBI_IDs_Jimmy_{}.txt'.format(fileNumber)
        
        print('Openeing {}'.format(fileName))

        substrate_data = open(fileName,'r')
        next(substrate_data)
        

	#un-pickle the dictionary and create a local copy
        
        substrates = pickle.load(open('substrate_dict.py','rb'))
        
        

	#print how many entries are in the dictionary currently

        print('Current Entries: {}'.format(len(substrates)))
        

	#initilize a dictionary to which values that are not present in the current version of the substrate data dictionary will be appended
	        
        new_entries = {}


	#go through the data file and check for any new proteins and their associated info

        for line in substrate_data:

                data = line.split('\t')

                if not substrates.has_key(data[0]):

                        new_entries[data[0]] = data[1:]

			#message stating the identifier of the protein that was added
            #print('{} added to dictionary.'.format(data[0]))


	#append the new entries to the substrate dictionary
        substrates.update(new_entries)


	# print the new number of entries and dump the dictionary back to the file
	print('\nUpdate Complete: {} total entries\n'.format(len(substrates)))
        pickle.dump(substrates,open('substrate_dict.py','wb'))



# A unit test written to check that the updates were in fact pickled

def test_dict():

	substrate = pickle.load(open('substrate_dict.py','rb'))

        print('Testing for Completion:{} entries in dictionary'.format(len(substrate)))



if __name__ == "__main__":

    for i in range(2,23):
        update_dict(i)
        test_dict()
