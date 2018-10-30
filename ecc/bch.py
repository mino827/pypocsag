# BCH Code implementation, shortened BCH for varying memory word width 32, 64, 128 ,256 , t = 2 
#Author: Colin 
#Date created: Nov 20 2013 
# Functionality: Encoder : systematic encoding, appending parity bits to the original information bits 
# all code vectors are interchangeably used with descending polynomial vector repr and number[bin/deci]. Two reprs: vector and number
# all polynomials are represented in descending order.  
# minimal polynomial table . Dependent on gfield.py module 

from ecc import gfield
import random 
import copy 
# minimal polynomial table. minpol[m][i] returns the number repr of minpoly for m_i(X) of GF(2^m). 
# used for syndrome calculation later 
minpol = {3:{1:11, 3:0b1011},
		4:{1:19,3:0b11111},
		5:{1:37,3:0b111101},
		6:{1:67, 3:0b1010111},
		7:{1:137, 3:0b10001111},
		8:{1:285, 3:0b101110111},
		9:{1:529, 3:0b1001011001}}

def CreateMessage(length):
	"create a information vector of length given, e.g. [0,0,1,1]"
	msg = []
	for i in range(0, length):
		letter = random.choice([0,1])
		msg.append(letter)

	return msg 

def encode(m,vin):
	"encode a input vector vin, and append r check bits , making a encoded word. r=2m, m determines GF(2^m)"
	#construct a FiniteField of order 2^m first 
	a = gfield.FiniteField(m)
	#minipol_1 = minpol[m][1] 
	#minipol_3 = minpol[m][3] 
	t = 2
	r = m*t # number of parity bits 
	# genpoly = a._poly_mul(minipol_1, minipol_3) # generator polynomial = LCM{m1,m3} 
	#genpoly = poly_mul(minipol_1, minipol_3) # generator polynomial = LCM{m1,m3} 
	genpoly = 0b11101101001
	# numerized_vin = a._vec2num(vin) # numerized_vin is a number 
	numerized_vin = gfield.vec2num(vin) # numerized_vin is a number 
	numerized_vin <<= r # X^r*u(x): left shift by r bits, i.e. adding trailing 0s for parity bits plugin 

	# rem = a._poly_div(numerized_vin, genpoly) 
	rem = gfield.poly_div(numerized_vin, genpoly) 
	# the number of leading 0 padding in the appendage parity vector 
	n_pad = r - len(bin(rem)) + 2 
	# trailing check bits vector . append it to original info vector
	v_checkbits = [0]*n_pad + list(a.showvector(rem))
	# append r parity bits to the original info vector 
	enc = v_checkbits 

	return enc 



def noise(vin, erate, nerr):
	"""corrupt the input vector vin by flipping at most nerr bits at a likelihood of erate , return a corrupted vector 
	CAUTION: vin is a list type(mutable), so it must be COPIED instead of assigned (which is python is 
	simply name-binding !!"""
	assert filter(lambda x:x in (0,1),vin)==vin  # ensure it's over GF(2) 
	noisy_v = copy.copy(vin) 
	count = 0 
	for i in range(len(noisy_v)):
		e = random.random()
		if e < erate:
			noisy_v[i] = int(not(noisy_v[i]))  # flip the bit 
			count += 1 

			if count == nerr:
				break 
		else:
			continue 

	return noisy_v 


def syndrome(m,vin):
	"""compute syndromes s1,s2,s3,s4 based on GF(2^m) operations 
	note s2 and s4 are implied with the knowledge of s1. since s4=s2^2, s2=s1^2 
	so only s1 and s3 are needed to compute error locator polynomial """
	a = gfield.FiniteField(m) 
	# numerize the input vector (list) to rec(received code as a number)
	if isinstance(vin, int):
		rec = vin 
	else:
		# rec = a._vec2num(vin) # convert to integer form 
		rec = gfield.vec2num(vin) # convert to integer form 
	# wrap_pol1 = a._poly_div(rec, minpol[m][1]) 
	wrap_pol1 = gfield.poly_div(rec, minpol[m][1]) 
	s1 = a.substitute(fin = wrap_pol1, powr=1) # S1 as a number 
	# wrap_pol3 = a._poly_div(rec, minpol[m][3]) 
	wrap_pol3 = gfield.poly_div(rec, minpol[m][3]) 
	s3 = a.substitute(fin = wrap_pol3, powr = 3) 


	return (s1,s3)  


def errorLocator(m,S1,S3):
	a = gfield.FiniteField(m) 
	# in case of s1=s3=0 : no error 
	if S1==0: # no error 
		return (0,0)
	A1 = S1 
	S1_cube = a.mul(S1, a.mul(S1,S1) ) 
	A2 = a.div(a.add(S3, S1_cube), S1)

	return (A1,A2)  


def chiensearch(m,A1,A2):
	"""chien search root in the GF(2^m) for error locator poly [A2,A1,1]
	input A1,A2,m ; all as number 
	output: roots=[root1,root2] so that alpha^root1 and alpha^root2 are roots 
	NOTE: must deal with 1 error case TODO """
	assert A1 # ensure at least 1 error exists 
	# assert A2 # ensure A1,A2 are non-zero ,i.e. errors exist 
	reg1 = A1 
	reg2 = A2 
	a = gfield.FiniteField(m)
	roots = []
	if a.add(reg1,reg2) == 1:
		roots.append(0) 
	for i in range(1,a.order - 1):
		reg1 = a.mul(reg1, 2)
		reg2 = a.mul(reg2, 4)
		summ = a.add(reg2, reg1) 
		if summ == 1: roots.append(i) 
	return roots 
		# CAUTION: deal with 1-error case 
	# if len(roots) <=2:
	# 	return roots 
	# else:
	# 	return (-1,-1) # decoding failure , more than 2 errors . this is source of decoding failure  

def errorPoly(m,A1,A2):
	"""given a error locator poly [A2,A1,1] in GF(2^m), find the error pattern poly e(X)
	a wrapper function of Chien search 
	return: a error polynomial as a number 
	NOTE: DEAL WITH 1-error case!! """
	if (A1,A2) == (0,0): 
		epoly = 0
		return epoly  # 0 error 
	a = gfield.FiniteField(m) 
	r_index = chiensearch(m, A1, A2) # r_index is a list of root indices 

	# def inverse(order,ind):
	# 	"get the inverse of a element index in GF(order), returns index "
	# 	inv = 0 if ind == 0 else order-1-ind 

	# 	return inv 

	inverse = lambda order,ind: 0 if ind == 0 else order-1-ind # inverse of a GF(order) element as index  




	if len(r_index) == 1: # 1 error 
		
		e_index = inverse(a.order, r_index[0])
		epoly = 1<<e_index 
		return epoly 
	elif len(r_index) == 2:  # 2 error 
		# inverse the roots index to get error location index 	
		e_index_1 = inverse(a.order, r_index[0])
		e_index_2 = inverse(a.order, r_index[1]) 
		epoly = (1<<e_index_1) ^ (1<<e_index_2)   
		return epoly  
	
	else:   # more than 3 err 
		epoly = -1 
		return epoly # decoding fails . -1 is a alarm status code more than 2 err


def correct(vin,errpol):
	"""
	input: vin(number), errpol (number) if not number but vector, convert to number first 
	output: corrected vector 
	"""
	assert isinstance(vin, list)  # ensure vin is a vector 
	assert errpol >= 0 # errpol must be non-negative. 
	if errpol == 0: # no error ! 
		return vin   
	

	n = len(vin) # n is code-length 
	
	rec = int(''.join(map(str, vin)), 2) # convert vin(list) to a number (int) 
	
	if isinstance(errpol, list):  #ensure to convert errpol to a number too
		err = int(''.join(map(str, errpol)), 2) 
	else:
		err = errpol 

	out_num = rec ^ err # c = r + e mod 2

	# convert out_num to vector form 

	n_pad = n - (len(bin(out_num)) - 2) 

	# vectorize out_num 
	out_vec = [0]*n_pad + list(map(int, bin(out_num)[2:]) ) 

	return out_vec 












