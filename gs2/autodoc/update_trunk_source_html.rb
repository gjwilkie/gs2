#!/home/edmundhighcock/Tools/bin/bin/ruby
require 'cgi'
require 'fileutils'
require 'pp'
source_dir = Dir.pwd + '/../temp'
html_dir = Dir.pwd + '/documentation/source'

`echo 'I was run' > null`
out_files = {}
function_references = {}
FileUtils.makedirs(source_dir)
FileUtils.makedirs(html_dir)	
FileUtils.makedirs(html_dir + '/geo')	
FileUtils.makedirs(html_dir + '/utils')	

Dir.chdir(source_dir) do
#	`svn checkout svn://gk.umd.edu/usr/local/repo/gs2/trunk/`
# 	Dir.chdir('trnk') do
		['.', 'utils', 'geo'].each do |dir| 
			Dir.chdir(dir) do
				Dir.entries(Dir.pwd).each do |file|	
					sub_dir = dir == "." ? "" : dir + "/"		
					out_file = "#{sub_dir}#{file}.html"
					(next) if not File.file? file or ['Makefile', 'README', 'test_os'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file))
					case File.extname(file)
					when "fpp"
						syntax = " -s=f90"
					else
						syntax = " "
					end
#					puts file
					out_files[html_dir + "/" + out_file] =  %x[highlight -H -a #{syntax} -i #{file}  --style lucretia --inline-css -K 12 -k Monaco -l].gsub(/(\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*subroutine\s*(?:\<\/span\>\s*)?(?:\<span[^>]*\>\s*)?)(\w+)(\s*.{40})/m){$2; function_references[$2] = "#{out_file}\##$2"; %[#$1<a name="#$2"></a>#$2#$3]}
				end
			end

		end
# 	end
end

File.open("function_references.rb", 'w'){|ref_file| ref_file.puts function_references.pretty_inspect} 

#puts function_references.pretty_inspect

#function_references.dup.each do |function, reference|
#	puts reference
#	function_references[function] = reference.gsub(/(?<double>[^\/]+\/)/){puts $~; gets; $~[:double]}
#end

#http://gyrokinetics.sourceforge.net/autodoc/
out_files.each do |out_file, text|
	puts out_file
	File.open(out_file, 'w'){|write_file| write_file.puts text.gsub(/(call(?:(?:(?:\s*\<\/span\>\s*)?(?:\s*\<span[^>]*\>\s*)?)|(?:\s+)))(\w+)/m){function_references[$2] ? %[#$1<a href= "#{Dir.pwd}/documentation/#{function_references[$2]}">#$2</a>] : "#$1#$2"} }
end

#puts "Content-type: text/plain; charset=iso-8859-1"

puts "Source Code html has been updated. Press your browser back button to continue"
