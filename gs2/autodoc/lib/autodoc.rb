#!/home/edmundhighcock/Tools/bin/bin/ruby

require 'cgi'
require 'fileutils'
require 'pp'
$not_found = []
class String
	attr_accessor :prefered_case
	def correct_case
		@prefered_case ||= :downcase
		self.send(@prefered_case.to_sym)
	end
	ADOCTOKENS = [['_ADOCDOT_', '.'], ['_ADOCFORWARDSLASH_', '/'], ['_ADOCPERCENT_', '%'], ['_ADOCHASH_', '#'], ['_ADOCCOLON_', ':'], ['_ADOCHYPHEN_', '-'], ['_ADOCSPACE_', ' ']]
	def from_adoc
		out = self.dup
		ADOCTOKENS.each{|pair| out.gsub!(pair[0], pair[1])}
		out
	end
	def to_adoc
		out = self.dup
		ADOCTOKENS.each{|pair| out.gsub!(pair[1], pair[0])}
		out
	end

end

# Autodoc is a class which generates documentation automatically from Fortran 95 source. Users write a small configuration file and run it using Ruby. As well as producing stylishly formatted HTML documentation Autodoc, using Highlight, also produces a fully hyperlinked HTML version of the source code which makes browsing the code extremely easy.
#
# Autodoc documents the structure of the code even without any user documentation. However, the user can document individual modules, functions and subroutines by placing comments surrounded by <tt><doc></doc></tt> braces. Comments placed at the beginning of the module before any subroutine definitions will document the module, and comments placed within the definition of a subroutine or function will document that subroutine or function.
#
# This documentation is aimed primarily at developers of the class. If you just want to use it have a look at the web page: 

class Autodoc
	class ParsingError < StandardError
	end
	
	SCRIPT_PATH = File.dirname(__FILE__) # :nodoc:
	

	FORTRAN_INTRINSIC = ["I", "ABORT", "ABS", "ACCESS", "ACHAR", "ACOS", "ACOSH", "ADJUSTL", "ADJUSTR", "AIMAG", "AINT", "ALARM", "ALL", "ALLOCATED", "AND", "ANINT", "ANY", "ASIN", "ASINH", "ASSOCIATED", "ATAN", "ATAN", "ATANH", "BESSEL", "BESSEL", "BESSEL", "BESSEL", "BESSEL", "BESSEL", "BIT", "BTEST", "C", "C", "C", "C", "C", "C", "CEILING", "CHAR", "CHDIR", "CHMOD", "CMPLX", "COMMAND", "COMPLEX", "CONJG", "COS", "COSH", "COUNT", "CPU", "CSHIFT", "CTIME", "DATE", "DBLE", "DCMPLX", "DFLOAT", "DIGITS", "DIM", "DOT", "DPROD", "DREAL", "DTIME", "EOSHIFT", "EPSILON", "ERF", "ERFC", "ERFC", "ETIME", "EXIT", "EXP", "EXPONENT", "FDATE", "FLOAT", "FGET", "FGETC", "FLOOR", "FLUSH", "FNUM", "FPUT", "FPUTC", "FRACTION", "FREE", "FSEEK", "FSTAT", "FTELL", "GAMMA", "GERROR", "GETARG", "GET", "GET", "GETCWD", "GETENV", "GET", "GETGID", "GETLOG", "GETPID", "GETUID", "GMTIME", "HOSTNM", "HUGE", "HYPOT", "IACHAR", "IAND", "IARGC", "IBCLR", "IBITS", "IBSET", "ICHAR", "IDATE", "IEOR", "IERRNO", "INDEX", "INT", "INT", "INT", "IOR", "IRAND", "IS", "IS", "ISATTY", "ISHFT", "ISHFTC", "ISNAN", "ITIME", "KILL", "KIND", "LBOUND", "LEADZ", "LEN", "LEN", "LGE", "LGT", "LINK", "LLE", "LLT", "LNBLNK", "LOC", "LOG", "LOG", "LOG", "LOGICAL", "LONG", "LSHIFT", "LSTAT", "LTIME", "MALLOC", "MATMUL", "MAX", "MAXEXPONENT", "MAXLOC", "MAXVAL", "MCLOCK", "MCLOCK", "MERGE", "MIN", "MINEXPONENT", "MINLOC", "MINVAL", "MOD", "MODULO", "MOVE", "MVBITS", "NEAREST", "NEW", "NINT", "NOT", "NULL", "OR", "PACK", "PERROR", "PRECISION", "PRESENT", "PRODUCT", "RADIX", "RAN", "RAND", "RANDOM", "RANDOM", "RANGE", "REAL", "RENAME", "REPEAT", "RESHAPE", "RRSPACING", "RSHIFT", "SCALE", "SCAN", "SECNDS", "SECOND", "SELECTED", "SELECTED", "SELECTED", "SET", "SHAPE", "SIGN", "SIGNAL", "SIN", "SINH", "SIZE", "SIZEOF", "SLEEP", "SNGL", "SPACING", "SPREAD", "SQRT", "SRAND", "STAT", "SUM", "SYMLNK", "SYSTEM", "SYSTEM", "TAN", "TANH", "TIME", "TIME", "TINY", "TRAILZ", "TRANSFER", "TRANSPOSE", "TRIM", "TTYNAM", "UBOUND", "UMASK", "UNLINK", "UNPACK", "VERIFY", "XOR"] # :nodoc:
	
	attr_accessor :code_website, :code_name
	
	# A description of the code which may contain HTML.
	
	attr_accessor :code_description
	
	#  A message which appears on the homepage of the code documentation. It can contain HTML tags.
	
	attr_accessor :welcome_message
	
	# Should Autodoc produce highlighted and hyperlinked HTML versions of the source code files? Default is true.
	
	attr_accessor :produce_highlighted_source_code
	
	attr_accessor :external_modules # :nodoc:
	attr_accessor :external_globals # :nodoc:
	attr_accessor :strict # :nodoc:

	# Any subdirectories in the source folder to include
	
	attr_accessor :subdirectories 

	
	# An array of filenames which Autodoc should not document.
	
	attr_accessor :ignore_files 
	
	# A hash of names and links which become tabs in the documentation. E.g. <tt>{"Wiki"=>"wiki-url"}</tt>
	
	attr_accessor :custom_tabs
	
	# A Proc object which takes the HTML output for each module page and performs any operation on it before returning the output. Allows the user to perform arbitrary operations on the documentation.
	
	attr_accessor :customize_documentation
	
	# A Proc object which takes the HTML output for each formatted source code file and performs any operation on it before returning the output. Allows the user to perform arbitrary operations on the highlighted source code.
	
	attr_accessor :customize_highlighting
	
	# If true, Autodoc will add hyperlinks to the GNU documentation for the Fortran intrinsic functions in the source code HTML. Makes Autodoc run more slowly.
	
	attr_accessor :document_fortran_intrinsic
	
	# A hash containing the hyperlinks to the functions in the source code. Structure is <tt> {[module_name, function_name] => "link"} </tt>
	
	attr_reader :function_references
	
	# A hash containing all the data that Autodoc has obtained by analysing the source code modules.
	
	attr_reader :modules
	
	# A hash containing information about which other modules and subroutines a particular subroutine or function uses.
	
	attr_reader :uses 
	
	# Create a new instance of Autodoc, which will analyse the source to be found in <tt>source_dir</tt> and generate documentation in <tt>html_dir</tt>
	
	def initialize(source_dir, html_dir)
		@source_dir, @html_dir = source_dir, html_dir
		@function_references = FORTRAN_INTRINSIC.inject({}) do |hash, func|
			hash[['fortran_intrinsic', func.correct_case]] = "http://gcc.gnu.org/onlinedocs/gcc-4.4.3/gfortran/#{func.upcase}.html##{func.upcase}"
			hash
		end #eval(File.read("function_references.rb"))
		@documentation_references = {}
		@out_pages = {}
		@highlighted_files = {}
		@subdirectories = [] #options[:sub_directories] ? options[:sub_directories] + ['.'] : ['.']
		@modules = {}
		@uses = {}
		@custom_tabs = {}
		@external_modules = []
		@external_globals = []
		@files_to_highlight = {}
		@ignore_files = []
		@strict = false
	end
	
	# Conditions under which a source file will not be documented.
	
	def exclude?(file, subdir_with_slash) #:doc:
		return (@ignore_files.include? "#{subdir_with_slash}#{file}" or not File.file? file or ['Makefile', 'README'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in', '.o', '.a', '.mod'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file)) or (false and not File.read(file)  =~ /\<wkdoc\>/))
	end
	
	private :exclude?
	
	# Highlight all the analysed source code files.
	
	def highlight_source 
		@highlighted_source = true
		analyse_source unless @analysed_source
		puts
# 		each_source_file do |file, subdir_with_slash|
		@files_to_highlight.each do |(file, subdir_with_slash), text|
			if @customize_highlighting
				highlighted =   @customize_highlighting.call(highlight_file(file, subdir_with_slash, text))
			else
				highlighted =  highlight_file(file, subdir_with_slash, text)
			end
			@highlighted_files["source/#{subdir_with_slash}#{file}.html"] = highlighted
		end
# 		FileUtils.rm(@source_dir + '/highlight.css')
		puts "\033[1A\033[KHighlighting: done"

# 		exit
	end
	
	# Highlight an individual source code file. First  Highlight is used to generate formatted HTML from the source code file. When the source code file was analysed any tokens due to be hyperlinked were replaced by a specially formatted token that Autodoc will recognise but which will not be altered by Highlight. These tokens are now replaced with the correct hyperlinks.
	# 
	# E.g.
	#	token_AUTODOC_NAME_name_AUTODOC_HREF_url_AUTODOC_DOCUMENTATION_documentation_url_AUTODOC_END
	#
	# becomes
	#
	# 	<a href="url">token</a><a name="name></a> (<a href="documentation_url">Documentation</a>)
	#
	# Note that all the URLs are escaped using special Autodoc syntax (see String#to_adoc).
	
	def highlight_file(file, subdir_with_slash, text) # :doc:
		puts "\033[1A\033[KHighlighting: #{file}"
		case File.extname(file)
		when "fpp"
			syntax = " -s=f90"
		else
			syntax = " "
		end
		FileUtils.makedirs('autodoc-temp')
		File.open('autodoc-temp/' + file, 'w'){|file| file.puts text}
		highlighted =  %x[highlight -H -a #{syntax} -i autodoc-temp/#{file}  --style kwrite --inline-css -K 12 -k Monaco -l]
		highlighted.gsub!(/(\w+?_AUTODOC_END)/) do
# 				p $1
			match = $1
			token = match.scan(/(^\w+?)_AUTODOC/).flatten[0]
			documentation = match.scan(/AUTODOC_DOCUMENTATION_(\w+?)_AUTODOC/).flatten[0]
			name = match.scan(/AUTODOC_NAME_(\w+?)_AUTODOC/).flatten[0]
# 				 if documentation
# 				(p 'start', match, token, name, documentation) #if documentation
			out = token
			href = match.scan(/AUTODOC_HREF_(\w+?)_AUTODOC/).flatten[0]
			if href 
				out = %[<a href="#{href.from_adoc}">#{out}</a>]
			end
			if name 
				out += %[<a name="#{name.from_adoc}"></a>]
			end
			
			if documentation
				out += %[ (<a href="#{documentation.from_adoc}">Documentation</a>)]
			end
# 			puts match, out
			out
		end

		highlighted
# 		p @function_references
# 		puts highlighted
# 		exit
# 		
# 		.gsub(/(call(?:(?:(?:\s*\<\/span\>\s*)?(?:\s*\<span[^>]*\>\s*)?)|(?:\s+)))(\w+)/m){@function_references[$2] ? %[#$1<a href= "#{@function_references[$2]}">#$2</a>] : "#$1#$2"} 
	end
	
	private :highlight_file
	
	# Hyperlink the given module. Within each subroutine or function replace all tokens which are actually other subroutines or functions used by the current one (as determined by the uses hash) with a hyperlink to that other subroutine or function. The hyperlinks are escaped using a special Autodoc notation (see Autodoc#highlight_file).
	
	def hyperlink_module(modname, modtext, file, subdir_with_slash) # :nodoc:
		sourcefolder = '../' * subdir_with_slash.scan(/\//).size
# 		puts "\n" * 8, modtext[0,20]
# 		p modname
		modname.gsub!(/_AUTODOC\w+_AUTODOC_END/i, '')
		modtext = each_subroutine_or_function(modtext) do |block, name|
			name.gsub!(/_AUTODOC\w+_AUTODOC_END/i, '')
			
# 			p modname, @modules[modname]
			#+ [['fortran_intrinsic', FORTRAN_INTRINSIC.map{|func| func.correct_case}]]
# 			p modname
			known_functions = ([[modname, (@modules[modname][:subroutines].keys or [])]] + (@uses[[modname, name]] or []).to_a)
			known_functions += [['fortran_intrinsic', FORTRAN_INTRINSIC.map{|func| func.correct_case}]] if @document_fortran_intrinsic
			
			known_functions.each do |usemod, funcs|
# 				p usemod
		                if funcs.size == 0
					funcs = @modules[usemod][:subroutines].keys if @modules[usemod]
		                end
# 				(p name, modname, known_functions, funcs; exit) if name =~ /betafun/ and modname =~ /geq/ and usemod =~ /splines/
		                funcs.each do |func|
# 					func = func.correct_case
					next if func == name
					if func =~ /\=\>/
						func, actual = func.split(/\s*\=\>\s*/)
					else
						actual = func
					end
# 					p func
					block.gsub!(Regexp.new("(\\W)#{func}(\\W)", Regexp::IGNORECASE)) do
						reference = @function_references[[usemod, actual]]
						beg = $1; en = $2
# 						(p name, usemod, func, actual, reference; exit) if name =~ /betafun/ and modname =~ /geq/ # unless reference or ["constants", "ifport", "hdf5", "mpi"].include? usemod #if func =~ /diameter/i
						if reference 
							unless reference =~ /^http/
								reference = sourcefolder + reference
							end
								%[#{beg}#{func}_AUTODOC_HREF_#{reference.to_adoc}_AUTODOC_END#{en}] 
								                             
						else
							%[#{beg}#{func}#{en}]
						end
					end
				end
			end
			block
		end[0]
		
	end
	
	private :hyperlink_module
	
	# Write the documentation to the @html_dir. First write all the hyperlinked source code files if hyperlinked source code has been switched on. Then write all the module documentation pages, and then main index page and the index of modules.
	
	def write_documentation
		FileUtils.makedirs(@html_dir)
		begin
			FileUtils.rm(@html_dir + '/styles.css')
			
		rescue
			
		end
		FileUtils.cp(SCRIPT_PATH + '/styles.css', @html_dir + '/styles.css')
		FileUtils.cp_r(SCRIPT_PATH + '/images', @html_dir + '/images')
		analyse_source unless @analysed_source
		highlight_source if @produce_highlighted_source_code and not @highlighted_source
		puts
		Dir.chdir(@html_dir) do 
			FileUtils.makedirs('source')
			Dir.chdir('source'){@subdirectories.each{|dir| FileUtils.makedirs dir}}
			@highlighted_files.each do |out_file, text|
				puts "\033[1A\033[KWriting source html: #{File.basename(out_file)}"
# 				puts out_file
				File.open(out_file, 'w'){|write_file| write_file.puts text}
			end
			puts "\033[1A\033[KWriting source html: done"
			puts
			@modules.each do |module_name, data|
				puts "\033[1A\033[KWriting doumentation: #{module_name}"
				File.open("#{module_name}.html", 'w')do |file| 
					mod = ModulePage.new(module_name, data, self).to_s
					mod = @customize_documentation.call(mod) if @customize_documentation 
					file.puts mod
				end
			end
			puts "\033[1A\033[KWriting doumentation: done"

			File.open("index.html", 'w')do |file| 
				index = IndexPage.new(self).to_s
				index = @customize_documentation.call(index) if @customize_documentation
				file.puts index
			end
			File.open("moduleindex.html", 'w')do |file| 
				index = ModuleIndexPage.new(self).to_s
				index = @customize_documentation.call(index) if @customize_documentation
				file.puts index
			end
			File.open("programindex.html", 'w')do |file| 
				index = ProgramIndexPage.new(self).to_s
				index = @customize_documentation.call(index) if @customize_documentation
				file.puts index
			end


		end
	end #def write_documentation

	# Iterate through each source file and pass each one to the block. subdir_with_slash is either <tt>""</tt> or <tt>subdir/</tt>, where subdir is a one of the subdirectories specified.
	
	def each_source_file(&block)
		Dir.chdir(@source_dir) do
		(@subdirectories + ['.']).each do |dir| 
			Dir.chdir(dir) do
				Dir.entries(Dir.pwd).each do |file|	
					subdir_with_slash = dir == "." ? "" : dir + "/"		
# 					out_file = "#{sub_dir}#{file}.html"
					(next) if exclude?(file, subdir_with_slash)
					yield(file, subdir_with_slash)
				end # Dir.entreis
			end #Dir.chdir
		end # @sub_directories.each
		end # Dir.chdir(source_dir)		
	end
	
	# Divide the text of the source file into sections, each section corresponding to one module. Each section starts after the word module, and ends after the <tt>end module name</tt> statement. Each section is passed to the block which may analyse or edit the block and but must return it. 
	#
	# Returns a list of <tt>[edited source, begining, ending]</tt> where beginning and ending are the sections of code not in any module definition (and which were consequently not passed to the block).
	
	def each_module(text, &block)
		module_blocks = (text).scan(/.+?(?:^\s*module(?![^\n]*procedure)|^\s*program|\Z)/im)
		beginning = module_blocks.shift; ending = module_blocks.pop
		all = beginning
		beginning =  yield(beginning, 'globals') # check for code not in modules
		ending =  yield(ending, 'globals') if ending
 		module_blocks.each do |modtext|
			modname = modtext.scan(/\A\s*(\w+)/).flatten[0].correct_case
			modname = 'program ' + modname if modtext =~ /^\s*end program/i
			(p file, modname, modtext; raise Autodoc::ParsingError) if modname =~ /procedure/i
			(puts modtext; raise ParsingError) if modname.length == 0
# 			p modname
			all += yield(modtext, modname)
		end
		[all+ending, beginning, ending] 
	end
	
	# Exactly the same as Autodoc#each_module except that this function separates the text of a module into individual sections each corresponding to the definition of a subroutine or a function

	def each_subroutine_or_function(text, &block)
		subroutine_blocks = (text).scan(/.+?(?:^\s*subroutine|^\s*(?:pure|integer)?\s*function|\Z)/im)
		beginning = subroutine_blocks.shift; ending = subroutine_blocks.pop
		all = beginning
		subroutine_blocks.each do |subtext|
			name = subtext.scan(/\A\s*(\w+)/)[0][0].correct_case
			all += yield(subtext, name)
		end
		[all+(ending or ""), beginning, ending] 
	end

	# Divide each source file into modules and pass each module to Autodoc#nalyse_module.
	
	def analyse_source
		@analysed_source = true
		puts
		each_source_file do |file, subdir_with_slash|
			puts "\033[1A\033[KAnalysing: #{subdir_with_slash}#{file}"

			begin
				text = File.read(file)
# 				next unless text  =~ /\<doc\>/
			rescue
				next
			end
			
			#(/.+?(?:\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*module(?:\s*\<\/span\>)?(?![^\n]*procedure)|\Z)/im)
			#/(.+^\s*module\s*(?!procedure)(\S+).*?\s+end\s+module)/im
			
			@files_to_highlight[[file, subdir_with_slash]] = (each_module("\n" + text + "\n\n") do |modtext, modname|
# 					pp modules[0]

				modtext = analyse_module(modname, modtext, file, subdir_with_slash)
# 				docfolder = '../' + '../' * subdir_with_slash.scan(/\//).size
				
			end)[0]

		end
		puts "\033[1A\033[KAnalysing: done"
		puts
		@files_to_highlight.each do |(file, subdir_with_slash), text|
			puts "\033[1A\033[KHyperlinking: #{subdir_with_slash}#{file}"

			@files_to_highlight[[file, subdir_with_slash]] = (each_module("\n" + text  + "\n\n") do |modtext, modname|
# 					pp modules[0]

				modtext = hyperlink_module(modname, modtext, file, subdir_with_slash)
# 				docfolder = '../' + '../' * subdir_with_slash.scan(/\//).size
			end)[0]
		end
				puts "\033[1A\033[KHyperlinking: done"
	
			
		File.open('uses.rb', 'w'){|file| file.puts @uses.pretty_inspect}
		File.open('fr.rb', 'w'){|file| file.puts @function_references.pretty_inspect}
		@modules['globals'][:description] = "Globally available functions and subroutines"

	end # analyse_source
	
	# Analyse the text of the module. Scan the module preamble for any documentation in closed in <tt><doc></doc></tt> braces. For each interface, subroutine or function the module defines add the correct reference to @function_references. For each subroutine and function, work out which other subroutine functions from other modules the current one uses and add that information to the variable @uses. For each subroutine and function, scanned for any documentation and add it to the <tt> @modules[module][:subroutines][subroutine]</tt> hash.
	
	
	def analyse_module(modname, modtext, file, subdir_with_slash) # :doc:


		@modules[modname] = {}
		@modules[modname][:file_name] = file
		@modules[modname][:subroutines] ||= {}
		@modules[modname][:subdir_with_slash] = subdir_with_slash unless modname == 'globals'
		docfolder = '../' + '../' * subdir_with_slash.scan(/\//).size
		documentation = docfolder + "#{modname}.html"
		modtext.sub!(/(\A\s*)(\w+)/){"#$1#$2_AUTODOC_NAME_#{modname.to_adoc}_AUTODOC_DOCUMENTATION_#{documentation.to_adoc}_AUTODOC_END"} unless modname == 'globals' 
		if @produce_highlighted_source_code
			interfaces = modtext.scan(/^\s*interface[ \t]+(\w+)/i).flatten
			(p modname, file, interfaces; raise ParsingError) if (["interface", "contains"] - interfaces).size < 2
			interfaces.each do |interface|
				@function_references[[modname, interface.correct_case]] = %[#{subdir_with_slash}#{file}.html##{modname}%#{interface.correct_case}] 
			end
		end
		modtext.gsub!(/(^\s*interface[\t ]+)(\w+)/i){iname = modname + '%' + $2; "#$1#$2_AUTODOC_NAME_#{iname.to_adoc}_AUTODOC_END"}
# 		puts modtext
# 		subroutine_blocks = modtext.scan(/.+?(?:^\s*subroutine(?![^\n]*procedure)|\Z)/im)
# 		beginning = subroutine_blocks.shift; ending = subroutine_blocks.pop
		
		modtext, beginning, ending =  each_subroutine_or_function(modtext) do |block, name|
			function_call = block.scan(/(\A.+(?:\&\s)?.+)/)[0][0].sub(/[\&\n]/, '')
			documentation = docfolder + "#{modname}.html##{name}"
			block.sub!(/(\A\s*)(\w+)/){"#$1#$2_AUTODOC_NAME_#{(modname +'%' + name).to_adoc}_AUTODOC_DOCUMENTATION_#{documentation.to_adoc}_AUTODOC_END"}

			comments = block.scan(/\<doc\>(.*?)\<\/doc\>/m).map{|match| match[0]}
			comments = comments.map do |comment|
				comment.gsub(/\n\s*\!/, '')
			end
# 			p 
			use_statements = block.scan(/(^\s*use[^&\n\r]+(?:(?:\&\s+)?[^&\n\r]+)+)/i)
# 			(p use_statements[0][0]; exit) if use_statements.find{|statements| statements[0] =~ /\&/} 
			@uses[[modname, name]] = {}
			use_statements.each do |statement_list|
# 				(p statement_list; exit) if  statement_list.grep(/constants.*nc_loop_fullmom/m)[0]
				words = statement_list[0].gsub(/\&/, '').scan(/\w+(?:\s*\=\>\s*)?\w+/).map{|w| w.correct_case}
# 				(p words; exit) if words.grep(/nc_loop_fullmom/)[0]
				words.delete("use"); words.delete("only")
				use_mod = words.shift
				use_subroutines = words
				@uses[[modname, name]][use_mod] ||= []
				@uses[[modname, name]][use_mod] += words
				@uses[[modname, name]][use_mod].uniq!
# 				p words; exit
			end
			@modules[modname][:subroutines][name] = [function_call, comments]
			@function_references[[modname, name]] = %[#{subdir_with_slash}#{file}.html##{modname}%#{name}] if @produce_highlighted_source_code
			@documentation_references[[modname, name]] = %[#{modname}.html##{name}]
			block
# 			beginning += block
		end
		@modules[modname][:description] = beginning.scan(/\<doc\>(.*?)\<\/doc\>/im).map{|match| match[0].gsub(/\n\s*\!/, '')}.join("\n")
# 		p beginning == modtext, ending
		return modtext #beginning + (ending or "")
 # 		exit if modname == "splines"
	end
	
	private :analyse_module
	
	# A class which knows how to generate the documentation HTML.
	
	class Page
		
		# Output the documentation HTML.
		
		def to_s
			<<EOF
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>#{@autodoccer.code_name} Documentation</title>
<meta name="keywords" content="" />
<meta name="description" content="" />
<link href="styles.css" rel="stylesheet" type="text/css" />
</head>
<body>
<div id="main"><!-- start header -->
<div id="header">
 <div id="logo">
	<h1>#{top_label}</h1>
	<h2><a href="http://www.metamorphozis.com/" id="metamorph">Designed by Metamorphosis Design,</a></h2><h2><a href="http://autodoc.rubyforge.org/" id="metamorph">produced by Autodoc</a></h2>
	</div>
	<div id="menu">
		<ul>
			<li><a href="index.html">Home</a></li>
			<li><a href="programindex.html">Programs</a></li>
			<li><a href="moduleindex.html">Modules</a></li>
			#{@autodoccer.custom_tabs.inject("") do |str, (name, ref)|
				str + %[<li><a href="#{ref}">#{name}</a></li>\n]
			end}
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
			#{program_sidebar}
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
		
		# The title in the top right-hand corner of the web page.
		
		def top_label
			%[<a href = "#{@autodoccer.code_website or "#"}">#{@autodoccer.code_name} Documentation</a>]
		end
		
		# A list of modules in the small div at the side of the web page.
		
		def module_sidebar
				keys = @autodoccer.modules.keys.sort.find_all{|mod| not mod =~ /program\s/i}
			<<EOF
				<h2>Module Index</h2>
				<ul class="back_title">
				#{keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="#{keys[0]}.html">#{keys[0]}</a></li>]) do |str, name|
					str + %[\n\t\t\t<li><a href="#{name}.html">#{name}</a></li>]
				end}
				</ul>
EOF
		end #def module_sidebar
		
		# A list of programs in the small div at the side of the web page.
		
		def program_sidebar
				keys = @autodoccer.modules.keys.find_all{|mod| mod =~ /program\s/i}
			<<EOF
				<h2>Program Index</h2>
				<ul class="back_title">
				#{keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="#{keys[0]}.html">#{keys[0]}</a></li>]) do |str, name|
					str + %[\n\t\t\t<li><a href="#{name}.html">#{name}</a></li>]
				end}
				</ul>
EOF
		end #def module_sidebar

		def subroutine_sidebar # :nodoc:
			return ""
		end #def subroutine_sidebar

	end # class page
	
	# A class which knows how to generate the HTML for the pages which document the individual modules.
	
	class ModulePage < Page
		def initialize(module_name, data, autodoccer)
			@module_name, @module_file_name, @subdir_with_slash, @subroutines, @autodoccer = module_name, data[:file_name], data[:subdir_with_slash], data[:subroutines], autodoccer
			@description = data[:description]
			 @function_references = autodoccer.function_references
		end
		
		# Generates the documentation for each subroutine or function in the main section of the module page. This consists of any comments made by the user, as well as information about which other functions this subroutine or function uses, a link to the source HTML and the call prototype.
		
		def subroutine_div(name, function_call, comments)
			lines = comments.map do |comment|
				line = "<li>#{comment}</li>"
				@function_references.each do |(modname, name), reference|
					next unless modname == name
					line.gsub(/name/, 
						%[<a href="#{reference}">name</a>])
				end
				line
			end
			uses = []
			use_string = ""
			if @autodoccer.uses[[@module_name, name]] and @autodoccer.uses[[@module_name, name]].size > 0
				@autodoccer.uses[[@module_name, name]].each do |usemod, words|
					use = %[ <li> <a href="#{usemod}.html">#{usemod}</a>, ]
					if words.size > 0
						use += "<small>only: </small>"  + words.inject("") do |str, word|
							usename, actual = word.split(/\s*\=\>\s*/)
							(actual = usename; usename = nil) unless actual
							str +=  %[name =&gt; ] if usename
							if @autodoccer.modules[usemod] and @autodoccer.modules[usemod][:subroutines].keys.include? actual
								str + %[ <a href="#{usemod}.html##{actual}">#{actual}</a>, ]
							else
								str + %[ #{actual}, ]
							end
						end
					else
						use += "<small>all</small>"
					end
					use += "</li>"
					uses.push use
				end
				use_string = %[<div class="notes">Uses:<ul>#{uses.join("\n")}</ul></div>]
			end
							
							
							
		return <<EOF
			#{@function_references[[@module_name, name]] ? %[<h2 class="title"><a href="source/#{@function_references[[@module_name, name]]}" name="#{name}">#{name}</a></h2>] : %[<h2 class="title">#{name}</h2>]} 
			
		<div class="notes"><small>Call Prototype:</small> #{function_call}</div>

			<div class="bullets"><!--<small>Comments:</small>--><ul>
				#{lines.join("\n\t\t\t")}
			</ul></div>
						#{use_string}
			 #{@function_references[[@module_name, name]] ? %[<div class="notes"><small><a href="source/#{@function_references[[@module_name, name]]}">View Source</a></small></div>] : %[]} 
			 <br>
EOF
		
		end #subroutine_div
		
		# Generate the content for the main section of the page. Start off with a list of all the functions and subroutines defined in the module and then has the documentation for individual subroutines (see ModulePage#subroutine_div).
		
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
	<!--<h2 class="title">Contains</h2>-->
	<div class = "entry">#{@subroutines.keys.sort.inject(""){|str, key| str += %[<a href="##{key}">#{key}</a>, ]; str}.chop.chop}</div><br>
#{@subroutines.inject("") do |str, (name, (function_call, comments))|
			str + subroutine_div(name, function_call, comments)
end}
	</div>
	<!-- end content -->
EOF
		end # def content
		
		# The information which appears in the blue box at the top of the page.
		
		def title_box
			<<EOF
<div id="box">
	<h1><a href="source/#@subdir_with_slash#{@module_file_name}.html##@module_name">#{@module_name =~ /program/i ? @module_name : "module " + @module_name} <small>(View Source)</small></a></h1>
<p>#{@description or "Brief description coming soon!"}</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s}</small></div>
</div>
EOF

# 	
	
		end #def title_box
# 		def subroutine_sidebar
# 				
# 				keys = @subroutines.keys.sort
# # 				(p @module_name, @subroutines; exit) if keys.size == 0
# 			<<EOF
# 				<h2>Subroutines</h2>
# 				<ul class="back_title">
# 				#{keys.size > 0 ? (keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="##{keys[0]}">#{keys[0]}</a></li>]) do |str, name|
# 					str + %[\n\t\t\t<li><a href="##{name}">#{name}</a></li>]
# 				end) : ""}
# 				</ul>
# EOF
# 		end #def subroutine_sidebar


	end # class ModulePage
	
	# This class generates the homepage for the documentation.
	
	class IndexPage < Page
		def initialize(autodoccer)
			@autodoccer =autodoccer
			@function_references = autodoccer.function_references
		end
		
		# Outputs the welcome message in the main section of the page.
		
		def welcome_message
			return <<EOF
			<h2 class="title">Welcome!</h2>
			<div class="entry">#{@autodoccer.welcome_message}</div>
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
		
		# Outputs the codename and the code description in the blue box at the top of the page. 
		
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">#{@autodoccer.code_name} Documentation Home Page</a></h1>
<p>#{@autodoccer.code_description}</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s}</small></div>
</div>
EOF
		end #def title_box


	end # class IndexPage

	# This class generates the page which has a list of all the modules and a short description of each module.
	
	class ModuleIndexPage < Page
		def initialize(autodoccer)
			@autodoccer =  autodoccer
			@modules = {}
			@autodoccer.modules.each do |name, data|
				@modules[name] = [data[:file_name], data[:description], data[:subdir_with_slash]]
			end
			 @function_references = autodoccer.function_references
		end
		
		def is_program_index? 
			false
		end
		
		private :is_program_index?
		
		# Generates a div which has the module name, short description and link to the source code.
		
		def module_div(name, file_name, description, subdir_with_slash)
# 			p file_name
		return <<EOF
<h2 class="title"><a href="#{name}.html" name="#{name}">#{name}</a></h2>
		<div class="entry"><small>Source File:</small><a href = "source/#{subdir_with_slash}#{file_name}.html"> #{file_name}</a></div>
		<div class="entry">#{description}</div><br>
EOF
		
		end #module_div
		
		# The content of the main section of the page.
		
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{@modules.inject("") do |str, (name, (file_name, description, subdir_with_slash))|
			if (name =~ /program/i and is_program_index?) or (not name =~ /program/i and not is_program_index?) 
				str + module_div(name, file_name, description, subdir_with_slash)
			else 
				str
			end	
			
end}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">Modules</a></h1>
<p>A list of modules with short descriptions.</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s} using Autodoc</small></div>
</div>
EOF
		end #def title_box
		def subroutine_sidebar
			return ""
		end #def subroutine_sidebar


	end # class ModuleIndexPage
	
	class ProgramIndexPage < ModuleIndexPage
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">Programs</a></h1>
<p>A list of available programs with short descriptions.</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s} using Autodoc</small></div>
</div>
EOF
		end
		def is_program_index? 
			true
		end
		
		private :is_program_index?

	end

end #class Autodoc
