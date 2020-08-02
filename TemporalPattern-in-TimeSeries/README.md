#### Background:
Time series is an ubiquitous type of data, appearing in many of the real-world problems, e.g. sensor data from vehicles. Identifying temporal patterns (e.g. subsequences)
in the time series to gain insight is an important part of time series analysis.

#### Problem definition:
Given a list of integers (time series data) in a text file separated by space.
1. Find the consecutive, non-empty subsequence with the highest _sum_. For example, , 3 -5 1 2 -1 4 -3 1 -2. The sum of all 9 integers is 0. 
=> the highest sum of two consecutive numbers is 3. For three consecutive numbers, it is 5. For a subsequence of four consecutive numbers is 6. 
=> So the highest sum of a subsequence is 6.
2. Instead of considering the whole sequence, find the highest possible sum from a subsequence, e.g. the highest sum of any subsequence of 10 integers from the sequence of 20 numbers.
3. Instead of direct _sum_ of the values we want to find the highest _sum_ of the absolute values of the differences of neighbouring pairs. For example,
absolute values of the differences of each pair would be `8 6 1 3 5 7 4 3` and therefore the output should be `16` for `n = 4` as `-1 4 -3 1` result in the absolute differences of `5 7 4` which adds up to `16`.

#### How to use:
All sequences are given with individual files placed in the 'data' folder. For problem 1 & 2, use the word 'values' in the argument when calling find_subsequence.py. On the other hand,
use the word 'differences' for problem-3. The command to run the program would be 'python find_subsequence.py inputFileLocation sequenceOrSubsequenceLength values_or_differences'.
For example, if you want get the maximum possible sum of a list given in 'input_1.txt', type in the command line "python find_subsequence.py data/input_1.txt".
Similarly, you can type "python find_subsequence.py data/input_1.txt 5 values" or "python find_subsequence.py data/input_1.txt 5 differences". Please have a look at few more example below:

\$ python find_subsequence.py data/input_1.txt 9 values
\>> 6
\$ python find_subsequence.py data/input_1.txt 4 differences
\>> 16
\$python find_subsequence.py data/input_2.txt 10 values
\>> 27
\$python find_subsequence.py data/input_2.txt 5 differences
\>> 58
\$python find_subsequence.py data/input_3.txt 30 values
\>> 44
\$python find_subsequence.py data/input_3.txt 10 differences
\>> 40
