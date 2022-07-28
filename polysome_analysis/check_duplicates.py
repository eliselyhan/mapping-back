with open("run_data.star") as f:
	lines = f.readlines()
	words = []
	for line in lines:
		if '5913-2_L2_ts002' in line:
			words.append(line)
			print(words)
			words_set = set(words)
print(len(words_set))

