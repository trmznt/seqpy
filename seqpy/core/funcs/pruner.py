
from seqpy import cerr, cexit

import numpy as np
import scipy.spatial.distance
import allel

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

	scoring_index = sorted( (s, i) for i, s in enumerate(score) )
	compress_index = [True] * len(scoring_index)

	# calculate r^2
	r = allel.rogers_huff_r(genotypes.to_n_alt())
	r_2 = scipy.spatial.distance.squareform( r**2 )

	# walk through index
	for i in range(N):
		idx = scoring_index[i]
		if not compress_index[i]:
			continue
		for _, j in scoring_index[i+1:]:
			if r_2[i,j] > threshold:
				compress_index[j] = False

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

	scoring_index = sorted( (s, i) for i, s in enumerate(score) )
	compress_index = [True] * len(scoring_index)

	# calculate r^2
	cerr('I: generating r^2 matrix')
	r = allel.rogers_huff_r(genotypes.to_n_alt())
	r_2 = scipy.spatial.distance.squareform( r**2 )

	# walk through index
	cerr('I: scanning r^2 matrix')
	for i in range(N):
		idx = scoring_index[i]
		if not compress_index[i]:
			continue
		for _, j in scoring_index[i+1:]:
			if r_2[i,j] > threshold:
				# check if this is a CDS region
				if positions[j][4] and positions[j][4] != positions[i][4]:
					continue
				compress_index[j] = False

	return compress_index
