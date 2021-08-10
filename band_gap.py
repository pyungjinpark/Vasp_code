#python for calculating band gap in insulator
def f():
	Positive_list = []
	Negative_list = []
	reading_file = open('band.txt','r')
	band = 0              
	band = reading_file.readlines()
	del band[0] #delete first line #kval total
	
	for line in band:
		line = line.strip()
		if len(line) != 0:
			line = line.replace('       ',' ')
			data = float(line.split()[1])
			if data > 0:
				Positive_list.append(data)
			elif data < 0:
				Negative_list.append(data)
	print('Check the band.dat file applied fermi level')

	print('max_negative:',max(Negative_list))
	print('min_positice:',min(Positive_list))
        Band_gap = min(Positive_list)-max(Negative_list)
        if Band_gap:
                print('band_gap:',format(Band_gap,'.4f'))
	else:	
		print('No band_gap')

        reading_file.close()


f()
