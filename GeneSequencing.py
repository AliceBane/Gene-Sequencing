# Alice Giola
# CS 4412
# Project 4 - Gene Sequencing
#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:
	def __init__(self):
		pass

	# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
	# handle to the GUI, so it can be updated as you find results, _banded_ is a boolean that tells
	# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
	# how many base pairs to use in computing the alignment
	def align(self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		results = []

		for i in range(len(sequences)):
			sequences[i] = sequences[i][:align_length] if len(sequences[i]) > align_length else sequences[i]

		for i in range(len(sequences)):
			jresults = []
			for j in range(len(sequences)):
				if (j < i):
					s = {}
				else:
					if banded:
						alignment1, alignment2, score = self.align_banded(sequences[i], sequences[j])
					else:
						alignment1, alignment2, score = self.align_unrestricted(sequences[i], sequences[j])
					s = {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
					table.item(i, j).setText('{}'.format(int(score) if score != math.inf else score))
					table.update()
				jresults.append(s)
			results.append(jresults)
		return results

	# The function align_unrestricted takes two strings s1 and s2 as input and initializes a 2D table (referred to
	# as table) with dimensions len(s2)+2 by len(s1)+2. The +2 is to allow for extra rows and columns at the top
	# and left sides of the table that will be used for initialization purposes later on.
	def align_unrestricted(self, s1, s2):
		table = Table(len(s2) + 2, len(s1) + 2)
		table.set_value(1, 1, TableValue(0))

		for char in range(len(s1)):
			value = TableValue((char + 1) * INDEL)
			table.set_value(1, char + 2, value)
			table.set_value(0, char + 2, s1[char])
		for char in range(len(s2)):
			value = TableValue((char + 1) * INDEL)
			table.set_value(char + 2, 1, value)
			table.set_value(char + 2, 0, s2[char])

		for y in range(2, len(s1) + 2):
			for x in range(2, len(s2) + 2):
				value_diagonal = table.get_value(x - 1, y - 1)
				value_left = table.get_value(x - 1, y)
				value_vert = table.get_value(x, y - 1)
				value_current = table.get_value(x, y)
				score_diagonal = value_diagonal.get_value()
				score_left = value_left.get_value() + INDEL
				score_vert = value_vert.get_value() + INDEL

				if table.get_value(x, 0) == table.get_value(0, y):
					score_diagonal += MATCH
				else:
					score_diagonal += SUB

				# Calculate the scores for each cell in the table using the three possible
				# moves (diagonal, left, and up).
				if score_left < score_vert and score_left < score_diagonal:
					value_current.set_value(score_left)
					value_current.set_pointer([x - 1, y])
				elif score_vert < score_diagonal and score_vert < score_left:
					value_current.set_value(score_vert)
					value_current.set_pointer([x, y - 1])
				elif score_diagonal < score_vert and score_diagonal < score_left:
					value_current.set_value(score_diagonal)
					value_current.set_pointer([x - 1, y - 1])

				# Handle ties
				elif (score_left == score_vert) or (score_left == score_diagonal):
					value_current.set_value(score_left)
					value_current.set_pointer([x - 1, y])
				elif score_vert == score_diagonal:
					value_current.set_value(score_vert)
					value_current.set_pointer([x, y - 1])
				else:
					value_current.set_value(score_diagonal)
					value_current.set_pointer([x - 1, y - 1])

		# Find the optimal alignment and score using the table.
		alignment1, alignment2 = self.find_alignment(table, len(s2) + 1, len(s1) + 1)
		score = table.get_value(len(s2) + 1, len(s1) + 1).get_value()

		# Return the optimal alignment and score.
		return alignment1, alignment2, score

	# The align_banded function takes in two strings s1 and s2 and returns the aligned sequences and
	# their alignment score. Huge thanks to Grant for trying to help me understand how to make it run in less than
	# 10 seconds and to my boyfriend for helping as well. I still couldn't figure out why some of the cells are off
	# by 4, but the final result is correct.
	def align_banded(self, s1, s2):
		band_distance = MAXINDELS
		if abs(len(s1) - len(s2)) > band_distance:
			return "", "", float('inf')  # Early exit if alignment is impossible within the band

		len_s1, len_s2 = len(s1), len(s2)
		band_width = 2 * band_distance + 1
		dp_table = [[float('inf')] * band_width for _ in range(len_s2 + 1)]
		backtrace = [[None] * band_width for _ in range(len_s2 + 1)]

		# Base case initialization for the DP table and backtrace paths
		dp_table[0][band_distance] = 0  # Starting point
		for i in range(band_distance + 1, band_width):
			dp_table[0][i] = dp_table[0][i - 1] + 5
			backtrace[0][i] = (0, i - 1)

		# Fill the DP table and backtrace within the band
		for i in range(1, len_s2 + 1):
			start = max(band_distance - i, 0)
			end = 6
			for j in range(start, end + 1):
				s1_index = i - band_distance + 2
				s2_index = i + j - band_distance - 1

				# Calculate costs for each possible move
				diag = float('inf')
				if 0 <= s1_index < len_s1 and 0 <= s2_index < len_s2:
					diag = dp_table[i - 1][j] + (MATCH if s1[s1_index] == s2[s2_index] else SUB)
				elif 0 <= s1_index < len_s1 and s2_index < 0:
					diag = dp_table[i - 1][j] + SUB
				left = dp_table[i][j - 1] + INDEL if j > 0 else float('inf')
				up = dp_table[i - 1][j + 1] + INDEL if j < end else float('inf')

				# Choose the best move
				min_cost = min(diag, left, up)
				dp_table[i][j] = min_cost
				if min_cost == diag:
					backtrace[i][j] = (i - 1, j)
				elif min_cost == left:
					backtrace[i][j] = (i, j - 1)
				else:  # up move
					backtrace[i][j] = (i - 1, j + 1)

		# Backtracking to reconstruct the alignment
		alignment1, alignment2 = '', ''
		i, j = len_s2, band_distance
		while i > 0 and j is not None:
			if backtrace[i][j] is None:
				break  # Safety check to avoid attempting to backtrack from an uninitialized cell
			prev_i, prev_j = backtrace[i][j]

			# Reconstruct alignment based on the move made to reach this cell
			s1_index = i + j - band_distance - 1
			if j == prev_j + 1 or j == prev_j:  # Diagonal or left move
				alignment1 = s1[s1_index] + alignment1 if s1_index < len_s1 else '-' + alignment1
				alignment2 = s2[i - 1] + alignment2
			else:  # Up move
				alignment1 = '-' + alignment1
				alignment2 = s2[i - 1] + alignment2 if prev_i < i else '-' + alignment2

			i, j = prev_i, prev_j

		if len_s1 < len_s2:
			score = dp_table[len_s2][band_distance + 1]
		else:
			score = dp_table[len_s2][band_distance]
		# score = dp_table[len_s2][band_distance]  # Fetch the final score from the DP table
		return alignment1, alignment2, score

	# takes in a table, and the final x and y coordinates of the table. The method retrieves the alignment strings by
	# tracing the path from the bottom right value of the table to the top left value, based on the pointers stored
	# in each value. Depending on the direction of the pointer, the method adds characters to either the first or
	# second alignment string, or adds a gap to one of them. Finally, the method returns the two alignment strings.
	def find_alignment(self, table, x, y):
		alignment1 = ''
		alignment2 = ''
		cost = table.get_value(x, y)

		while cost.get_pointer() is not None:
			pointer_x, pointer_y = cost.get_pointer()
			# diagonal move
			if pointer_x == x - 1 and pointer_y == y - 1:
				alignment1 = table.get_value(0, y) + alignment1
				alignment2 = table.get_value(x, 0) + alignment2
			# vertical move
			elif pointer_x == x - 1:
				alignment1 = '-' + alignment1
				alignment2 = table.get_value(x, 0) + alignment2
			# horizontal move
			else:
				alignment1 = table.get_value(0, y) + alignment1
				alignment2 = '-' + alignment2

			x, y = pointer_x, pointer_y
			cost = table.get_value(x, y)

		return alignment1, alignment2


# A defined table value with a pointer location and a value inside of it.
# Used in conjunction with the value table array matrix above.
class TableValue:
	def __init__(self, starting_value, input_pointer=None):
		self.pointer = input_pointer
		self.value = starting_value

	def set_pointer(self, input_pointer):
		self.pointer = input_pointer

	def get_pointer(self):
		return self.pointer

	def set_value(self, input_value):
		self.value = input_value

	def get_value(self):
		return self.value


# Table class with methods for setting and getting values in a two-dimensional table, as well as a str method for
# printing the table. The set_value method sets the value at the given (x, y) coordinate in the table, while the
# get_value method returns the value at that coordinate.
class Table:
	def __init__(self, x_size, y_size):
		self.table = []
		for i in range(y_size):
			x_array = []
			for j in range(x_size):
				x_array.append(TableValue(""))
			self.table.append(x_array)

	# might need rework (not sure if direct input of values to the
	# table will work 100% percent of the time).
	def set_value(self, x_input, y_input, input_value):
		table = self.table
		table[y_input][x_input] = input_value

	def get_value(self, x_input, y_input):
		table = self.table
		return table[y_input][x_input]

	def __str__(self):
		temp_table = self.table
		output_string = ''

		for i in range(len(temp_table)):
			for j in range(len(temp_table[i])):
				output_string += str(temp_table[i][j]) + ' '
			output_string += "\n"
		return output_string
