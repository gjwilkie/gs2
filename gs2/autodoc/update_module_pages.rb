#!/home/edmundhighcock/Tools/bin/bin/ruby

require 'cgi'
require 'fileutils'
require 'pp'


class Autodoc
	def initialize(source_dir, html_dir, options={})
		@source_dir, @html_dir = source_dir, html_dir
		@function_references = eval(File.read("function_references.rb"))
		@out_pages = {}
		@sub_directories = options[:sub_directories] ? options[:sub_directories] + ['.'] : ['.']
	end
	def analyse_source
		Dir.chdir(source_dir) do
		@sub_directories.each do |dir| 
			Dir.chdir(dir) do
				Dir.entries(Dir.pwd).each do |file|	
					sub_dir = dir == "." ? "" : dir + "/"		
					out_file = "#{sub_dir}#{file}.html"
					(next) if not File.file? file or ['Makefile', 'README', 'test_os'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file)) 
					begin
						text = File.read(file)
						next unless text  =~ /\<wkdoc\>/
					rescue
						next
					end
					@subroutines = {}
					subroutine_blocks = text.split(/^\s*subroutine\s*/)
					subroutine_blocks.slice(1..(subroutine_blocks.size)).each do |block|
						name = block.scan(/(\A\w+)/)[0][0]
						p block
						puts "\n\n\n"
						function_call = block.scan(/(\A.+(?:\&\s)?.+)/)[0][0].sub(/[\&\n]/, '')
# 						function_call = block.scan(/(\A.*?(?:\(.*[\n]?.*\))?)/)[0][0].sub(/[\&\n]/, '')
						comments = block.scan(/\<wkdoc\>(.*?)\<\/wkdoc\>/m).map{|match| match[0]}
						comments = comments.map do |comment|
							comment.gsub(/\n\s*\!/, '')
						end
						@subroutines[name] = [function_call, comments]
					end
				end # Dir.entreis
			end #Dir.chdir
		end # @sub_directories.each
		end # Dir.chdir(source_dir)
	end # analyse_source
	def subroutine_div(name, function_call, comments)
		lines = comments.map do |comment|
			line = "<li>#{comment}</li>"
			function_references.each do |name, reference|
				line.gsub(/name/, 
					%[<a href="#{reference}">name</a>])
			end
			line
		end
		<<EOF
		<div class="entry"><small>Call Prototype:</small> #{function_call}</div>
			<ul>
				#{lines.join("\n\t\t\t")}
			</ul>
EOF
		
					end
end
# source_dir = ENV['PWD'] + '/../gs2_source'
source_dir = Dir.pwd + '/../temp'



`echo 'I was run' > null`
out_pages = {}
function_references = eval(File.read("function_references.rb"))


Dir.chdir(source_dir) do
#	`svn checkout svn://gk.umd.edu/usr/local/repo/gs2/trunk/`
# 	Dir.chdir('trunk') do
		['.', 'utils', 'geo'].each do |dir| 
			Dir.chdir(dir) do
				Dir.entries(Dir.pwd).each do |file|	
					sub_dir = dir == "." ? "" : dir + "/"		
					out_file = "#{sub_dir}#{file}.html"
					(next) if not File.file? file or ['Makefile', 'README', 'test_os'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file)) 
					begin
						text = File.read(file)
						next unless text  =~ /\<wkdoc\>/
					rescue
						next
					end
					subroutines = {}
					subroutine_blocks = text.split(/^\s*subroutine\s*/)
					subroutine_blocks.slice(1..(subroutine_blocks.size)).each do |block|
						name = block.scan(/(\A\w+)/)[0][0]
						p block
						puts "\n\n\n"
						function_call = block.scan(/(\A.+(?:\&\s)?.+)/)[0][0].sub(/[\&\n]/, '')
# 						function_call = block.scan(/(\A.*?(?:\(.*[\n]?.*\))?)/)[0][0].sub(/[\&\n]/, '')
						comments = block.scan(/\<wkdoc\>(.*?)\<\/wkdoc\>/m).map{|match| match[0]}
						lines = comments.map do |comment|
							comment = comment.gsub(/\n\s*\!/, '')
							comment = "<li>#{comment}</li>"
							function_references.each do |name, reference|
								comment.gsub(/name/, 
									%[<a href="#{reference}">name</a>])
							end
							comment
						end
						subroutines[name] = %[<div class="entry"><small>Call Prototype:</small> #{function_call}</div><ul>] + lines.join("\n") + "</ul>"
					end
					keys = subroutines.keys #.sort do |name1, name2|
# 						name1 <=> name2
# 					end
#.sub(/\..*$/, '')
					out_pages[file + ".html"] = keys.inject(<<EOF
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Metamorphosis Design Free Css Templates</title>
<meta name="keywords" content="" />
<meta name="description" content="" />
<link href="styles.css" rel="stylesheet" type="text/css" />
</head>
<body>
<div id="main"><!-- start header -->
<div id="header">
 <div id="logo">
	<h1><a href="#">metamorph_darkside</a></h1>
	<h2><a href="http://www.metamorphozis.com/" id="metamorph">Design by Metamorphosis Design</a></h2>
	</div>
	<div id="menu">
		<ul>
			<li><a href="#">Home</a></li>
			<li><a href="#">Blogs</a></li>
			<li><a href="#">Photos</a></li>
			<li><a href="#">About</a></li>
			<li><a href="#">Contact</a></li>
		</ul>
	</div>	
<!-- end header -->
</div>
<!-- start page -->
<div id="top"></div>	
<div id="page">
<div id="box">
	<h1><a href="#">Documentation for GS2 Module: #{file.sub(/\..*$/, '')}</a></h1>
<p>Brief description coming soon!</p>
</div>
	  <!-- start content -->
	<div id="content">
EOF
					                                        ) do |str, name|

						str + <<EOF
			#{function_references[name] ? %[<h2 class="title"><a href="source/#{function_references[name]}">#{name}</a></h2>] : %[<h2 class="title">#{name}</h2>]} 
			<div class="entry"><small>Last updated #{Time.now.to_s} using <a href="#">Autodoc</a></small></div>
		#{subroutines[name]}


EOF
# ===#{function_references[name] ? "[http://gs2wiki.edmundhighcock.com/trunk_source/#{function_references[name]} #{name}]" : name}===
# 
# #{subroutines[name]}

					end + <<EOF
	</div>
	<!-- end content -->
	<!-- start sidebar two -->
	<div id="sidebar2" class="sidebar">
			<h2>Categories</h2>
				<ul class="back_title">
					<li class="top"><a href="#">Aliquam libero</a></li>
					<li><a href="#">Consectetuer elit</a></li>
					<li><a href="#">Metus pellentesque</a></li>
					<li><a href="#">Suspendisse mauris</a></li>
				</ul>
				<h2>Archives</h2>
				<ul class="back_title">
					<li class="top"><a href="#">September</a> (23)</li>
					<li><a href="#">August</a> (31)</li>
					<li><a href="#">July</a> (31)</li>
					<li><a href="#">June</a> (30)</li>
					<li><a href="#">May</a> (31)</li>
				</ul>
				<h2>Lorem ipsum dolor </h2>
				<ul class="back_title">
					<li class="top"><a href="#">Metus pellentesque</a></li>
					<li><a href="#">Suspendisse mauris</a></li>
					<li><a href="#">Urnanet molestie semper</a></li>
					<li><a href="#">Proin orci porttitor</a></li>
				</ul>
	</div>
	<!-- end sidebar two -->
<div style="clear: both;">&nbsp;</div>

</div>
	<div id="bottom"></div>
<!-- end page -->
<!-- start footer -->
<div id="footer">
	 <p>Copyright &copy; 2009. <a href="#">Privacy Policy</a> | <a href="#">Terms of Use</a> | <a href="http://validator.w3.org/check/referer" title="This page validates as XHTML 1.0 Transitional"><abbr title="eXtensible HyperText Markup Language">XHTML</abbr></a> | <a href="http://jigsaw.w3.org/css-validator/check/referer" title="This page validates as CSS"><abbr title="Cascading Style Sheets">CSS</abbr></a></p> 
	<p>Design by <a href="http://www.metamorphozis.com/" title="Free Flash Templates">Free Flash Templates</a>
		</p>
</div>
</div>
<!-- end footer -->
</body>
</html>

EOF

					
# 					out_pages[file + ".html"] = keys.inject("==Member Subroutines==\n\n''Do not edit this section - it is generated automatically from the source code''.\n\n") do |str, name|
# 
# 						str + <<EOF
# ===#{function_references[name] ? "[http://gs2wiki.edmundhighcock.com/trunk_source/#{function_references[name]} #{name}]" : name}===
# 
# #{subroutines[name]}
# 
# EOF
# 					end + "''Do not add anything after this point - it will be deleted''\n\n[[Category:Modules]]" 
# 
				end
			end

		end
# 	end
end

Dir.chdir('documentation') do 
	  out_pages.each do |file_name, doctext|
		File.open(file_name, 'w'){|file| file.puts doctext}
	  end
end
puts out_pages.values

# Dir.chdir('/home/edmundhighcock/wikibot/pywikipedia/') do 
# 	out_pages.each do |page, doctext|
# 		%x[python2.5 replace.py -page:#{page} -regex "==Member Subroutines==.*(\n.*)+" "#{doctext.gsub(/\\/, '\\\\').gsub(/"/, '\\"')}" -always]
# 	end
# end    
# 
# #puts "Content-type: text/plain; charset=iso-8859-1"
# 
# puts "Source Code html has been updated. Press your browser back button to continue"
