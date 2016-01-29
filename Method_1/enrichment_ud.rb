require 'fileutils'

def median(arr)
	len = arr.length
	sorted = arr.sort
	median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
	return median
end

def median_line(lines)
	min_FDR = 1.5
	med_fc = 0
	aux = Array.new	
	gene = lines.first.downcase.split(/\t/).first
	if lines.length == 2
		lines.each{|line|
                        data = line.chomp.downcase.split(/\t/)
			if data[2].to_f < min_FDR
				med_fc = data[1].to_f
				min_FDR = data[2].to_f			
			end
		}
	else
		lines.each{|line|
			data = line.chomp.downcase.split(/\t/)
			aux << data[1].to_f
			min_FDR = [min_FDR,data[2].to_f].min
		}
		med_fc = median(aux)
	end
	return "#{gene}\t#{med_fc}\t#{min_FDR}"
end


def get_csv(file)
	f = File.open(file)
	f.gets
	genesets = Array.new
        gene_targets = Hash.new
	while (line = f.gets)
		data = line.downcase.chomp.split(/;/)
		genesets << data[0]
		gene_targets[data[2]] ||= Array.new
		gene_targets[data[2]] << data[0]
	end
	f.close
	return [genesets.uniq, gene_targets]

end


def get_gmx(file)
	f = File.open(file) 
	genesets = f.gets.chomp.split(/\t/)
	gene_targets = Hash.new

	while (line = f.gets)
		data = line.downcase.chomp.split(/\t/)
		i = 0
		data.each{|g|
			if !g.nil? and !g.empty?
				gene = g.strip	
				gene_targets[gene] ||= Array.new
				gene_targets[gene] << genesets[i]
			end
			i+=1
		}
	end
	f.close
	return [genesets,gene_targets]
end

def get_data(file)

	i = File.open(file)
	fl = i.gets.chomp
	content = Hash.new
	while (line = i.gets)
		data = line.downcase.split(/\t/)
		content[data[0]] ||= Array.new
		content[data[0]] << line.downcase.chomp
	end
	i.close
	return [fl,content]
end

if ARGV[0].nil? or ARGV[1].nil? or ARGV[2].nil? or ARGV[3].nil?
		puts "\n\tPerforms the miRNA enrichment in experimentally up- and donwregulated genes\n\n"
        puts "\tUsage: ruby enrichment_ud.rb <input_file> <pvalue_threshold> <down_threshold> <up_threshold>\n\n"
        puts "Instructions of use\n------------------\n"
		puts "\n1- A database of predicted interactions is required with the following format:\n\n"
		puts "\tmiRNA;gene;gene\n\t<miRNA 1>;<gene 1>;<gene 1>\n\t.\n\t.\n\t.\n\t<miRNA n>;<gene n>;<gene n>"
		puts "Its name must be 'union.csv' or changed in line number 99. See example given with the same name\n\n"
		puts "2- The <input file> must be a tab-delimited three-column file containing for every gene of interest its name,\n"
		puts "fold-change and p-value. The p-value can be set to 0 if those values are not available in the expermiental data.\n"
		puts "See example file named example.txt\n\n"
		puts "3- <p-value threshold>: higher limit of pvalue for a miRNA to be printed in the output\n\n"
		puts "4- <down_threshold>: higher limit for downregulation\n\n"
		puts "5- <up_threshold>: lower limit for upregulation\n\n"
		puts "6- The outputs are two files named:\n\n"
		puts "\t a) <input name>_results.txt:\n\n"
		puts "A tab-delimited file consisting of 19 columns:\n"
		puts "'miRNA': microRNA name\n'totalExpected': total number of genes in the predictive database used\n'expectedUp': number of upregulated genes in the input data\n"
		puts "'expectedDn': number of downregulated genes in the input data\n'totalObserved': total number of genes in the predictive database used\n'observedUp': number of genes in the database annotated as upregulated in the input file\n"
		puts "'observedDn': number of genes in the database annotated as downregulated in the input file\n'UpReg expected ratio': expectedUp\/totalExpected\n'DnReg expected ratio: expectedDn\/totalExpected'\n"
		puts "'UpReg observed ratio': observedUp\/totalObserved\n'DnReg observed ratio': observedDn\/totalObserved\n"
		puts "'UpReg confidence interval 95% lower'\n'UpReg confidence interval 95% upper'\n"
		puts "'Upreg pvalue': Wilson's p-value for enrichment among upregulated targets\n"
		puts "'DnReg confidence interval 95% lower'\n'DnReg confidence interval 95% upper'\n"
		puts "'Dnreg pvalue': Wilson's p-value for enrichment among downregulated targets\n"
		puts "'Permutation test Up pValue': p-value by a permutation test, for enrichment among upregulated targets\n"
		puts "'Permutation test Dn pValue': p-value by a permutation test, for enrichment among downregulated targets\n\n"
		puts "\t b) <input name>_genesets.txt:\n\n"
		puts "A tab-delimited file summarizing the targets predicted for every miRNA\n"
		puts "See example output files given\n\n"
        exit
end

# GMX format
#(genesets, gene_targets) = get_gmx("gene_set.gmx")

# miRNAs from predictive DB format
(genesets, gene_targets) = get_csv("union.csv")
(fl, content) = get_data(ARGV[0])
name = ARGV[0].split(/\./).first

o = File.open("#{ARGV[0].gsub(/\.txt/,"")}_genesets.txt","w")
o.write("#{fl}\t#{genesets.join("\t")}\n")

content.each_pair{|gene,lines|

	if lines.length > 1
		o.write(median_line(lines))
	else
		o.write(lines.first)
	end
	genesets.each{|m|
		if !gene_targets[gene].nil? and gene_targets[gene].include?(m)
			o.write("\t1") 
		else
			o.write("\t0")
		end
	}
	o.write("\n")
}
o.close


# Matlab invoke

mfile = File.open("script.m","w")
mfile.write("scriptcossgsea('#{name}_genesets.txt', '#{name}_results.txt', #{ARGV[1].to_f}, #{ARGV[2].to_f}, #{ARGV[3].to_f}, #{genesets.length + 3})")
mfile.close

`matlab -nojvm -nodesktop -nodisplay < script.m`

FileUtils.rm "script.m"
