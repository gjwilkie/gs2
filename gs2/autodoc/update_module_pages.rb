#!/home/edmundhighcock/Tools/bin/bin/ruby

require 'cgi'
require 'fileutils'
require 'pp'


class Autodoc
	attr_accessor :code_website, :code_name, :code_description, :welcome_message
	attr_reader :function_references, :modules
	def initialize(source_dir, html_dir, options={})
		@source_dir, @html_dir = source_dir, html_dir
		@function_references = eval(File.read("function_references.rb"))
		@out_pages = {}
		@highlighted_files = {}
		@sub_directories = options[:sub_directories] ? options[:sub_directories] + ['.'] : ['.']
		@title_prefix = (options[:title_prefix] or "Autodoc Documentation")
		@modules = {}
	end
	def exclude?(file)
		return (not File.file? file or ['Makefile', 'README', 'test_os'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file)))
	end
	def analyse_and_highlight_source_code
		Dir.chdir(@source_dir) do
		@sub_directories.each do |dir| 
			Dir.chdir(dir) do
				Dir.entries(Dir.pwd).each do |file|	
					sub_dir = dir == "." ? "" : dir + "/"		
					out_file = "#{sub_dir}#{file}.html"
					(next) if exclude? file
#					puts file
					@highlighted_files[@html_dir + "/source/" + out_file] =  analyse_and_highlight_file(file, out_file)
				end
			end

		end
		@highlighted_files.each do |out_file, text|
		puts out_file
			File.open(out_file, 'w'){|write_file| write_file.puts text.gsub(/(call(?:(?:(?:\s*\<\/span\>\s*)?(?:\s*\<span[^>]*\>\s*)?)|(?:\s+)))(\w+)/m){@function_references[$2] ? %[#$1<a href= "#@html_dir/source/#{@function_references[$2]}">#$2</a>] : "#$1#$2"} }
		end

	end
# end

# File.open("function_references.rb", 'w'){|ref_file| ref_file.puts function_references.pretty_inspect} 

#puts function_references.pretty_inspect

#function_references.dup.each do |function, reference|
#	puts reference
#	function_references[function] = reference.gsub(/(?<double>[^\/]+\/)/){puts $~; gets; $~[:double]}
#end

#http://gyrokinetics.sourceforge.net/autodoc/

	end #analyse_and_highlight_source_code
	def analyse_and_highlight_file(file, out_file) 
		case File.extname(file)
		when "fpp"
			syntax = " -s=f90"
		else
			syntax = " "
		end
		return %x[highlight -H -a #{syntax} -i #{file}  --style lucretia --inline-css -K 12 -k Monaco -l].gsub(/(\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*subroutine\s*(?:\<\/span\>\s*)?(?:\<span[^>]*\>\s*)?)(\w+)(\s*.{40})/m){$2; @function_references[$2] = "#{out_file}\##$2"; %[#$1<a name="#$2"></a>#$2#$3]}
	end
	def write_documentation
		analyse_documentation
		Dir.chdir(@html_dir) do 
			@modules.each do |module_name, data|
				File.open("#{module_name}.html", 'w')do |file| 
					mod = ModulePage.new(module_name, data, self)
					file.puts mod
				end
			end
			File.open("index.html", 'w')do |file| 
				index = IndexPage.new(self)
				file.puts index
			end
			File.open("moduleindex.html", 'w')do |file| 
				index = ModuleIndexPage.new(self)
				file.puts index
			end


		end
	end #def write_documentation

	def analyse_documentation
		Dir.chdir(@source_dir) do
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
					modules = text.scan(/(^\s*module\s*(\S+).*?\s*end\s+module)/m)
# 					pp modules[0]
					modules.each do |modtext, name|
						analyse_module_documentation(name, file, modtext)
					end

				end # Dir.entreis
			end #Dir.chdir
		end # @sub_directories.each
		end # Dir.chdir(source_dir)
	end # analyse_documentation
	def analyse_module_documentation(modname, file, modtext)
		@modules[modname] = {}
		@modules[modname][:file_name] = file
		@modules[modname][:subroutines] = {}
		subroutine_blocks = modtext.split(/^\s*subroutine\s*/)
		@modules[modname][:description] = subroutine_blocks[0].scan(/\<wkdoc\>(.*?)\<\/wkdoc\>/m).map{|match| match[0].gsub(/\n\s*\!/, '')}.join("\n")
		subroutine_blocks.slice(1..(subroutine_blocks.size)).each do |block|
			name = block.scan(/(\A\w+)/)[0][0]
# 						p block
# 			puts "\n\n\n"
			function_call = block.scan(/(\A.+(?:\&\s)?.+)/)[0][0].sub(/[\&\n]/, '')
# 						function_call = block.scan(/(\A.*?(?:\(.*[\n]?.*\))?)/)[0][0].sub(/[\&\n]/, '')
			comments = block.scan(/\<wkdoc\>(.*?)\<\/wkdoc\>/m).map{|match| match[0]}
			comments = comments.map do |comment|
				comment.gsub(/\n\s*\!/, '')
			end
			@modules[modname][:subroutines][name] = [function_call, comments]
		end
	end
	class Page
		def to_s
			<<EOF
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
	<h1>#{top_label}</h1>
	<h2><a href="http://www.metamorphozis.com/" id="metamorph">Design by Metamorphosis Design</a></h2>
	</div>
	<div id="menu">
		<ul>
			<li><a href="index.html">Home</a></li>
			<li><a href="moduleindex.html">Modules</a></li>
			<!--<li><a href="https://sourceforge.net/apps/mediawiki/gyrokinetics">Wiki</a></li>
			<li><a href="#">About</a></li>
			<li><a href="#">Contact</a></li>-->
		</ul>
	</div>	
<!-- end header -->
</div>
<!-- start page -->
<div id="top"></div>	
<div id="page">
#{title_box}
#{content}
	<!-- start sidebar -->
	<div id="sidebar2" class="sidebar">
			#{subroutine_sidebar}
			#{module_sidebar}
	</div>
	<!-- end sidebar -->
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
		end # def page
		def top_label
			%[<a href = "#{@autodoccer.code_website or "#"}">#{@autodoccer.code_name} Documentation</a>]
		end
		def module_sidebar
				keys = @autodoccer.modules.keys
			<<EOF
				<h2>Modules</h2>
				<ul class="back_title">
				#{keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="#{keys[0]}.html">#{keys[0]}</a></li>]) do |str, name|
					str + %[\n\t\t\t<li><a href="#{name}.html">#{name}</a></li>]
				end}
				</ul>
EOF
		end #def module_sidebar
		def subroutine_sidebar
			return ""
		end #def subroutine_sidebar

	end # class page
	class ModulePage < Page
		def initialize(module_name, data, autodoccer)
			@module_name, @module_file_name, @subroutines, @autodoccer = module_name, data[:file_name], data[:subroutines], autodoccer
			@description = data[:description]
			 @function_references = autodoccer.function_references
		end
		def subroutine_div(name, function_call, comments)
		lines = comments.map do |comment|
			line = "<li>#{comment}</li>"
			@function_references.each do |name, reference|
				line.gsub(/name/, 
					%[<a href="#{reference}">name</a>])
			end
			line
		end
		return <<EOF
			#{@function_references[name] ? %[<h2 class="title"><a href="source/#{@function_references[name]}" name="#{name}">#{name}</a></h2>] : %[<h2 class="title">#{name}</h2>]} 
			
		<div class="entry"><small>Call Prototype:</small> #{function_call} #{@function_references[name] ? %[<small><a href="source/#{@function_references[name]}">View Source</a></small>] : %[]} </div>
			<ul>
				#{lines.join("\n\t\t\t")}
			</ul><br>
EOF
		
		end #subroutine_div
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{@subroutines.inject("") do |str, (name, (function_call, comments))|
			str + subroutine_div(name, function_call, comments)
end}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">Module: #@module_name</a></h1>
<p>#{@description or "Brief description coming soon!"}</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s} using Autodoc</small></div>
</div>
EOF
		end #def title_box
		def subroutine_sidebar
				keys = @subroutines.keys
			<<EOF
				<h2>Subroutines</h2>
				<ul class="back_title">
				#{keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="##{keys[0]}">#{keys[0]}</a></li>]) do |str, name|
					str + %[\n\t\t\t<li><a href="##{name}">#{name}</a></li>]
				end}
				</ul>
EOF
		end #def subroutine_sidebar


	end # class ModulePage
	
	class IndexPage < Page
		def initialize(autodoccer)
			@autodoccer =autodoccer
			@function_references = autodoccer.function_references
		end
		def welcome_message
			return <<EOF
			<h2 class="title">Welcome!</h2>
			<p>#{@autodoccer.welcome_message}</p>
			<div class="entry"><small> This documentation has been automatically generated from the #{@autodoccer.code_name} source code by Autodoc</small></div>
EOF
		end # def welcome_message
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{welcome_message}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">#{@autodoccer.code_name} Documentation Home Page</a></h1>
<p>#{@autodoccer.code_description}</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s} using Autodoc</small></div>
</div>
EOF
		end #def title_box


	end # class IndexPage

	class ModuleIndexPage < Page
		def initialize(autodoccer)
			@autodoccer =  autodoccer
			@modules = {}
			@autodoccer.modules.each do |name, data|
				@modules[name] = [data[:file_name], data[:description]]
			end
			 @function_references = autodoccer.function_references
		end
		def module_div(name, file_name, description)
# 			p file_name
		return <<EOF
<h2 class="title"><a href="#{name}.html" name="#{name}">#{name}</a></h2>
		<div class="entry"><small>Source File:</small><a href = "source/#{file_name}.html"> #{file_name}</a></div>
		<div class="entry"><p>#{description}</p></div>
EOF
		
		end #module_div
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{@modules.inject("") do |str, (name, (file_name, description))|
			str + module_div(name, file_name, description)
end}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">Module: #@module_name</a></h1>
<p>A list of modules with short descriptions.</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s} using Autodoc</small></div>
</div>
EOF
		end #def title_box
		def subroutine_sidebar
			return ""
		end #def subroutine_sidebar


	end # class ModuleIndexPage

end #class Autodoc

autodoccer = Autodoc.new(Dir.pwd + '/../temp', Dir.pwd + '/documentation', {:sub_directories => ['utils', 'geo']})
autodoccer.code_name = "GS2"
autodoccer.code_website = "http://gyrokinetics.sourceforge.net"
autodoccer.code_description = "GS2 is a gyrokinetic flux tube initial value turbulence code which can be used for fusion or astrophysical plasmas."
autodoccer.analyse_and_highlight_source_code
autodoccer.write_documentation