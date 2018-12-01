
from seqpy import cerr, cexit

import numpy as np
import scipy.spatial.distance
import allel

def arrange_index( score ):
	""" return index based on descending order of score """
	scoring_index = sorted( (s, i) for i, s in enumerate(score), reverse=True )
	index = np.array( idx for _, idx in scoring_index, dtype=np.int32 )
	return index

def calculate_r_2( genotypes ):
	cerr('I: calculating Rogers-Huff r^2...')
	r = allel.rogers_huff_r(genotypes.to_n_alt())
	r_2 = scipy.spatial.distance.squareform( r**2 )
	return r_2


def prune_1(genotypes, threshold=0.5, score=None):
	""" prune by r^2 with score as priority,
		returning indexing array
	"""

	N = len(genotypes)

	if score == None:
		# we use MAC as default score
		score = np.min( genotypes.count_alleles(), axis=1)

	if N != len(score):
		cexit('E: length of genotypes !+ length of score!')

	index = arrange_index( score )
	compress_index = np.ones( len(index), dtype=np.int8 )

	# calculate r^2
	r_2 = calculate_r_2( genotypes )

	# walk through index
	for i in range(N):
		if not compress_index[i]:
			continue
		for j in index[i:]:
			if r_2[i,j] > threshold:
				compress_index[j] = 0

	return compress_index


def prune_2(genotypes, positions, threshold=0.5, score=None):
	""" prune by r^2 except on CDS, only CDS within the same segment/region will be pruned
	"""

	if score == None:
		# we use MAC as default score
		score = np.min( genotypes.count_alleles(), axis=1)

	N = len(genotypes)

	if N != len(score) or N != len(positions):
		cexit('E: length of genotypes != length of score nor positions!')

	index = arrange_index( score )
	compress_index = np.ones( len(index), dtype=np.int8 )

	# calculate r^2
	r_2 = calculate_r_2( genotypes )

	# walk through index
	for i in range(N):
		if not compress_index[i]:
			continue
		for j in scoring_index[i:]:
			if r_2[i,j] > threshold:
				# check if this is a CDS region
				if positions[j][4] and positions[j][4] != positions[i][4]:
					continue
				compress_index[j] = 0

	return compress_index


def prune_21(genoarray, positions, threshold=0.5):
	""" prune by r^2 except on CDS, only CDS within the same segment/region will be pruned,
		with priority given to highest non-missingness - maf
	"""

	ca = genoarray.count_alleles()
	total_ca = np.sum( ca )
	maf = np.min( ca / total_ca )
	completeness = 1 - ca.count(-1)/total_ca
	score = np.array( completeness, maf )

	return prune_2(genoarray, positions, threshold, score)
