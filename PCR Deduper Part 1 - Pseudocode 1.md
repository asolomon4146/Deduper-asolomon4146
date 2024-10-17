## Definition of problem
PCR duplicate definition: An artifact in DNA sequencing where identical DNA fragments are over-represented due to excessive amplification during PCR. It is important to remove these PCR duplicates because they introduce bias which can affect conclusions.
## Requirements to be a PCR dupe
Same UMI
Same alignment position (after adjusting for soft clipping)
Same chromosome
Same Strand
# Pseudocode
```
Initialize empty dict "dupe_dict"
initialize an empty string named "key"

Open input file and output file simultaneously, use "with open" method
	for line in enumerate(input file): # reads input sam file line by line
		If starts with @, write to out.sam # writes the headers
		else:
			split line by tab to make each line into a tab delimited list of columns # line = line.split()
			extract last part of column 1 to extract UMI.
			if UMI is in known list of UMIs, add it to the string named "key"
			extract strand information by checking bitwise flag of bit 16 and append to "key" separated by a semicolon
			extract CHROM from column 3 and append to "key" separated by a semicolon
			Check if soft clipped by checking if bit 16 is true in column 5 CIGAR
				if ((flag & 16) == 16) :
					rev_comp = True
					This makes 
				Use regex to extract the number which precedes the S in the CIGAR string or until it reaches another character or whitespace.
					Append it to the string "key"
			if string named "key" does not exist in dupe_dict:
				add it to dupe_dict as a key and add the line as the value.
	for i in dupe_dict.values(): 
		iterate through the list of dupe_dict values and write each one to the out.sam file.
```
