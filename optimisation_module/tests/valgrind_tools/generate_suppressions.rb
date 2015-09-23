#p ARGV;exit
system("valgrind --gen-suppressions=all --error-limit=no #{ARGV.join(' ')} 2> .supptmpfile")
text = File.read(".supptmpfile")
text = text.gsub(/^==.*$/, '').scan(/\{.*?\}/m).uniq.join("\n")
File.open('valgrind.supp', 'w'){|f| f.puts text}

