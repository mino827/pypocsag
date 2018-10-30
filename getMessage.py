#!/usr/bin/env python3

#Author zezadas
#Credits nunohumberto

import sys
from ecc import bch, gfield
import random
NRZ_I = True

VERBOSE = False
INJECT_ERROR = False
BCH_M = 5

codeword_size=32 #32bits
frame_size=2*codeword_size #1frame=2words
batch_size=codeword_size + (8*frame_size) #sc + 8 frames
iso7bit_size=7 #ISO 7-bits ISO 646
num_bit_size=4

preamble = "10"*288 #576bits 1 and 0
sync_word = "01111100110100100001010111011000" #0x7CD215D8
idle_word = "01111010100010011100000110010111" #0x7A89C197

addr=""
msg_words=""
msg_type=""

def print_verbose(msg):
    if VERBOSE:
        print(msg)

def loadFile(filename):
	f = open(filename, 'rb')
	return f.read().strip()


def process_batch(payload):
	global msg_words
	global addr
	global msg_type
	sc= payload[:32]
	if sc != sync_word:
		print("incorrect sync word")
		sys.exit()
	for word_idx in range(1,int(len(payload)/codeword_size)):
		word = payload[codeword_size*word_idx:codeword_size*(word_idx+1)]
		if word[0] == "0": ##adr
			if word == idle_word:
				print_verbose("-----------------")
				print_verbose("New idle word: " + word)
			else:
				addr = int(word[1:19]+format(int((word_idx-1)/2), '03b'),2)
				msg_type= word[19:21]
				print_verbose("-----------------")
				print_verbose("New address word: " + word)
		if word[0] == "1":
			print_verbose("-----------------")
			print_verbose("New message word: " + word)
			checkMessage(word)
			msg_words+=word[1:21]


def getMessage(bits):

	raw_data = bits.decode('utf-8')

	if NRZ_I:
		raw_data = ''.join([["0","1"][i=="0"] for i in raw_data])

	print_verbose(raw_data)

	if raw_data[0:576] != preamble:
		print("incorrect preamble")
		sys.exit()

	payload = raw_data[576:]

	if len(payload) < batch_size:
		print("incorrect payload size %s" % len(payload))
		sys.exit()

	for batch_idx in range(int(len(payload)/batch_size)):
		print_verbose("new batch %s" % batch_idx)
		batch = payload[batch_size*batch_idx:batch_size*(batch_idx+1)]
		process_batch(batch)

	if msg_type=="11":
		msg=""
		for i in range(int(len(msg_words)/iso7bit_size)):
			msg+= chr(int(msg_words[iso7bit_size*i:iso7bit_size*(i+1)][::-1],2))
	elif msg_type=="00":
		msg=""
		dc=["<Spare>","<Urgent>","<space>","'","]","["]
		for i in range((len(msg_words)/num_bit_size)):
			num = int(msg_words[num_bit_size*i:num_bit_size*(i+1)][::-1],2)
			if num > 9:
				dc_idx= num-10
				num=dc[dc_idx]
			msg += str(num)
	else:
		print("msg type not recognized")
		sys.exit()
	print_verbose("-----------------")
	print("addr: %s - msg: %s" % (addr,msg))


# This method checks for errors in a received message.
def checkMessage(word):
	word_parity = int(word[-1]) # Get parity bit
	word_data_and_ecc = word[:-1] # Get data + BCH bits

	num_array = [int(a) for a in word_data_and_ecc] # Convert to number array

	int_data = num_array[0:21] # Obtain data array
	int_bch = num_array[21:31] # Obtain BCH array

    #for test purposes
	if INJECT_ERROR:	# Inject error in data array at a random position
		original = int_data[:]		# Keep original data to evaluate error correction
		generateSingleError(int_data)

	one_counter = sum(int_data+int_bch)		# Count number of ones to verify parity
	
	if one_counter & 1 != word_parity:
		print("Parity check: Failed.")
	else:
		print_verbose("Parity check: OK.")

	if VERBOSE:		# Show differences between expected BCH bits and received BCH bits
		expected_bits = ''.join([str(c) for c in bch.encode(BCH_M, int_data)])
		obtained_bits = ''.join([str(c) for c in int_bch])
		print_verbose("Expected BCH bits: " + expected_bits)
		print_verbose("Received BCH bits: " + obtained_bits)


	(s1,s3) = bch.syndrome(BCH_M, int_data+int_bch) 

	if s1 == 0:		# Check if errors are present
			print_verbose('No error.')
			return

	# Error correction		
	# ----------------
	(A1,A2) = bch.errorLocator(BCH_M, s1, s3) 
	err_pol = bch.errorPoly(BCH_M, A1, A2) 

	if err_pol < 0:
		print_verbose('Error decoding.')
		

	corr_v = bch.correct(int_data+int_bch, err_pol) 
	# ----------------

	if INJECT_ERROR:
		print_verbose('Original  data: ' + str(original))

	print_verbose('Corrupted data: ' + str(int_data))
	print_verbose('Recovered data: ' + str(corr_v[0:21]))



# This method generates a single error at a random position in an array of bits
def generateSingleError(bits):
	index = int(round(random.random()*(len(bits)-1)))
	bits[index] = 0 if bits[index] else 1
	return bits 


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Usage: " + sys.argv[0] + " <input bit sequence file>")
		sys.exit(1)
	bits = loadFile(sys.argv[1])
	getMessage(bits)
