package Parsers::YYLexer;
#
# $RCSfile: YYLexer.pm,v $
# $Date: 2017/01/12 18:59:59 $
# $Revision: 1.12 $
#
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2017 Manish Sud. All rights reserved.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#

use strict;
use Carp;
use Exporter;
use Scalar::Util ();
use Parsers::Lexer;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Parsers::Lexer Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyYYLexer';

# Class constructor...
sub new {
  my($Class, $Input,  @TokensSpec) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new($Input,  @TokensSpec);
  bless $This, ref($Class) || $Class;
  $This->_InitializeYYLexer();

  return $This;
}

# Initialize object data...
#
sub _InitializeYYLexer {
  my($This) = @_;

  # File generated containing mapping of token labels to token numbers by
  # running command byacc with -d option on a parser definition file.
  #
  # For example, "byacc -l -P -d -b Parser Parser.yy" would generate file
  # Parser.tab.ph, which might contain the following tokem name and values
  # for a parser for a simple calculator:
  #
  #  $NUMBER=257;
  #  $LETTER=258;
  #
  #
  $This->{YYTabFile} = undef;
  $This->{YYTabFilePath} = undef;

  # Mapping of token lables to token numbers...
  %{$This->{YYTabDataMap}} = ();

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...

  $ClassName = __PACKAGE__;
}

# Process tokens in YYTab file and load mapping of token labels to integers
# for return during YYLex method invocation...
#
# Notes:
#   . YYTabFile must be a complete path or available through @INC path in the
#     same directory where this package is located.
#   . Name of YYTabFile might start with any valid sub directory name in @INC
#     For example, "Parsers/<YYTablFile>" implies the tab file in parsers sub directory
#     under MayaChemTools lib directory as it would be already in @INC path.
#   . YYTabFile must be explicitly set by the caller. The default YYTabFile name,
#     y.tab.ph, generated by byacc is not used implicitly to avoid confusion among
#     multiple distinct instances of YYLexer.
#   . YYTabFile is generated by byacc during its usage with -d options and contains
#     mapping of token codes to token names/labels. YYLexer used this file to map
#     token lables to token codes before retuning token code and value pair back
#     to yyparse function used by byacc.
#   . User defined token numbers start from 257
#
#     The format of YYTabFile generted by byacc during generation of parser code in
#     Perl code is:
#
#     ... ...
#     $NUMBER=257;
#     $ADDOP=258;
#     $SUBOP=259;
#     ... ..
#
sub SetupYYTabFile {
  my($This, $YYTabFile) = @_;
  my($YYTabFilePath, $Line, $TokenLabel, $TokenNumber);

  $This->{YYTabFile} = undef;
  $This->{YYTabFilePath} = undef;
  %{$This->{YYTabDataMap}} = ();

  if (!defined $YYTabFile) {
    croak "Error: ${ClassName}->SetupYYTabFile: YYTabFile must be specified...";
  }
  $This->{YYTabFile} = $YYTabFile;

  if (-e $YYTabFile) {
    $YYTabFilePath = $YYTabFile;
  }
  else {
    ($YYTabFilePath) = grep {-f "$_/$YYTabFile"}  @INC;
    if (!$YYTabFilePath) {
      carp "Warning: ${ClassName}->SetupYYTabFile: YYTabFile, $YYTabFile,  can't be located in \@INC path: @INC...";
      return $This;
    }
    $YYTabFilePath = "${YYTabFilePath}/$YYTabFile";
  }

  $This->{YYTabFilePath} = $YYTabFilePath;

  open YYTABFILE, "$YYTabFilePath" or die "Couldn't open $YYTabFilePath: $_\n";
  while ($Line = <YYTABFILE>) {
    ($TokenLabel, $TokenNumber) = ($Line =~ /^\$(.*?)=(.*?);$/);
    if (!(defined($TokenLabel) && defined($TokenNumber))) {
      croak "Error: ${ClassName}->SetupYYTabFile: Couldn't extract token label and number from YYTabFile $YYTabFile at line: $Line...";
    }
    if (exists $This->{YYTabDataMap}{$TokenLabel}) {
      carp "Warning: ${ClassName}->SetupYYTabFile: Token lable, $TokenLabel, already defined in YYTabFile $YYTabFile...";
    }
    $This->{YYTabDataMap}{$TokenLabel} = $TokenNumber;
  }
  close YYTABFILE;

  return $This;
}

# Get next available token number and any matched text from input stream
# by either removing it from the input stream or simply peeking ahead.
#
# Supported mode values: Peek, Next. Default: Next
#
# Notes:
#   . Token label and value pairs returned by lexer, which can't be mapped to token
#     labels specified in YYTabFile are ignored.
#   . Token text of length 1 returned by lexer without a corresponding explicit token label,
#     which can't be mapped to a token number using Perl ord function, is ignored.
#
sub YYLex {
  my($This, $Mode) = @_;
  my($LexerToken, $TokenLabel, $TokenNumber, $TokenText);

  ($TokenLabel, $TokenNumber, $TokenText) = (undef) x 3;

  TOKEN: while (defined($LexerToken = $This->Lex($Mode))) {
    if (ref $LexerToken) {
      ($TokenLabel, $TokenText) = @{$LexerToken};
      if (exists $This->{YYTabDataMap}{$TokenLabel}) {
	$TokenNumber = $This->{YYTabDataMap}{$TokenLabel};
      }
      elsif ($TokenLabel =~ /^EOI$/i) {
	$TokenNumber = 0;
      }
    }
    else {
      $TokenText = $LexerToken;
    }

    # Check for any literals (+, - , = etc.) to generte token numbers...
    #
    if (!defined $TokenNumber) {
      if (length($TokenText) == 1 && ord $TokenText) {
	$TokenNumber = ord $TokenText;
      }
    }

    # Generate error message for no mapping to token numbers...
    if (defined $TokenNumber) {
      last TOKEN;
    }
    else {
      if (defined $TokenLabel) {
	carp "Warning: ${ClassName}->YYLex: Igorning token label, $TokenLabel, with matched text, $TokenText, returned by lexer and retrieving next available token or text. Token label couldn't be mapped to token numbers specified in YYTabFile generated from a parser defintion file using byacc. After updating parser definition file, a new YYTabFile containing entry for token label must be generated...";
      }
      else {
	carp "Warning: ${ClassName}->YYLex: Igorning token text, $TokenText, returned by lexer and retrieving next available token or text. Token text returned by lexer couldn't be mapped to token number using Perl ord function. After updating lexer token specifications and parser definition file, a new YYTabFile containing entry for a new token label to match unrecognized text must be generated...  ";
      }
      next TOKEN;
    }
  }

  if (!defined $LexerToken) {
    # Chained lexer returns undefined at end of input. So it's equivalent to EOI
    # token.
    if (exists $This->{YYTabDataMap}{EOI}) {
      $TokenLabel = "EOI";
      $TokenNumber = $This->{YYTabDataMap}{$TokenLabel};
      $TokenText = "0";
    }
    else {
      ($TokenLabel, $TokenNumber, $TokenText) = ("EOI", 0, "0");
    }
  }

  return ($TokenNumber, $TokenText);
}

# Get next available token number and text pair from input stream by removing it
# from the input stream...
#
sub Next {
  my($This) = @_;

  return $This->YYLex();
}

# Get next available token number and text pair from input stream by by simply
# peeking ahead and without removing it from the input stream...
#
sub Peek {
  my($This) = @_;

  return $This->YYLex('Peek')
}

# Return a curried verson of lexer: yyparse in parser generated by byacc expects it
# to call without passing any argument for the YYLexer object...
#
sub GetYYLex {
  my($This) = @_;

  return sub { my($Mode) = @_; $This->YYLex($Mode); };
}

# Is it a lexer object?
sub _IsYYLexer {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Return a string containing information about lexer...
sub StringifyYYLexer {
  my($This) = @_;
  my($YYLexerString);

  $YYLexerString = "YYLexer: PackageName: $ClassName; " . $This->_GetYYLexerInfoString();

  return $YYLexerString;
}

# Stringigy YYTabFile token name and value information...
#
sub _GetYYLexerInfoString {
  my($This) = @_;
  my($YYLexerInfoString, $TokenValue, $YYTabFile, $YYTabFilePath, $YYTabDataMapString);

  $YYTabFile = defined $This->{YYTabFile} ? $This->{YYTabFile} : 'None';
  $YYTabFilePath = defined $This->{YYTabFilePath} ? $This->{YYTabFilePath} : 'None';

  $YYLexerInfoString = "YYTabFile: $YYTabFile; YYTabFilePath: $YYTabFilePath";

  $YYTabDataMapString = "YYTabDataMap: None";
  if (keys %{$This->{YYTabDataMap}}) {
    my($TokenLabel, $TokenNumber);

    $YYTabDataMapString = "YYTabDataMap:";
    for $TokenLabel (sort keys %{$This->{YYTabDataMap}}) {
      $TokenValue = $This->{YYTabDataMap}{$TokenLabel};
      $YYTabDataMapString .= " ${TokenLabel}=${TokenValue}";
    }
  }

  $YYLexerInfoString .= "; $YYTabDataMapString; " . $This->_GetLexerInfoString();

  return $YYLexerInfoString;
}

1;

__END__

=head1 NAME

Parsers::YYLexer

=head1 SYNOPSIS

use Parseres::YYLexer;

use Parsers::YYLexer qw(:all);

=head1 DESCRIPTION

B<YYLexer> class provides the following methods:

new, GetYYLex, Next, Peek, SetupYYTabFile, StringifyYYLexer, YYLex

B<Parsers::YYLexer> class is derived from B<Parsers::Lexer> base class, which provides all
the underlying lexer functionality. B<YYLexer> class is designed to be used with
B<yyparse> code generated by running B<byacc> on a parsers defined using
parser definition B<ParserName.yy> file.

I<YYTabFile> containing mapping of token labels to integers must be explicitly
specified by the caller. This file is processed during new method invocation and
mapping of token labels to integers is loaded in a hash to be used later by B<YYLex>
method to return token number and text pairs to the parser.

=head2 METHODS

=over 4

=item B<new>

    $YYLexer = new Parsers::YYLexer($Input,  @YYLexerTokensSpec);

Using specified I<Input> and I<YYLexerTokensSpec>, B<new> method generates a new
B<YYLexer> and returns a reference to newly created B<YYLexer> object.

Examples:

    # Tokens specifications supplied by the caller. It's an array containing references
    # to arrays with each containing TokenLabel and TokenMatchRegex pair along with
    # an option reference to code to be executed after a matched.
    #
    @LexerTokensSpec = (
        [ 'LETTER', qr/[a-zA-Z]/ ],
        [ 'NUMBER', qr/\d+/ ],
        [ 'SPACE', qr/[ ]*/,
            sub { my($This, $TokenLabel, $MatchedText) = @_; return ''; }
        ],
        [ 'NEWLINE', qr/(?:\r\n|\r|\n)/,
            sub { my($This, $TokenLabel, $MatchedText) = @_;  return "\n"; }
        ],
        [ 'CHAR', qr/./ ]
    );

    # Input string...
    $InputText = 'y = 3 + 4';

    $YLexer = new Parsers::YYLexer($InputText,  @LexerTokensSpec);

    # Setup default token table file...
    $YYTabFile = "Parsers/SimpleCalcParser.tab.ph";
    $This->SetupYYTabFile($YYTabFile);

    # Process input stream...
    ($TokenNumber, $TokenText) = $YYLexer->Lex();
    print "TokenNumber: $TokenNumber; TokenText: $TokenText\n";

    # Input file...
    $InputFile = "Input.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $Lexer = new Parsers::YYLexer(\*INPUTFILE, @LexerTokensSpec);

    # Input file iterator...
    $InputFile = "TestSimpleCalcParser.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $InputIterator = sub { return <INPUTFILE>; };
    $Lexer = new Parsers::YYLexer($InputIterator, @LexerTokensSpec);

    # Usage with code generated by byacc from a parser definition
    # file SimpleCalcParser.yy...

    $InputText = "3 + 4 +6\nx=3\ny=5\nx+y\nx+z\n";

    $YYLexer = new Parsers::YYLexer($InputText,@LexerTokensSpec);

    $YYLex = $YYLexer->GetYYLex();

    $YYTabFile = "Parsers/SimpleCalcParser.tab.ph";
    $YYLexer->SetupYYTabFile($YYTabFile);

    $Debug = 0;
    $SimpleCalcParser = new Parsers::SimpleCalcParser($YYLex,
                            \&Parsers::SimpleCalcParser::yyerror, $Debug);

    $Value = $SimpleCalcParser->yyparse();
    print "Value = " . (defined($Value) ? "$Value" : "Undefined") . "\n";

=item B<GetYYLex>

    $YYLex = $YYLexer->GetYYLex();

Returns a curried verson of YYLexer as B<YYLex>: yyparse in parser generated by
byacc expects it to call without passing any argument for the I<YYLexer> object.

=item B<Next>

    ($TokenNumber, $TokenText) = $YYLexer->Next();

Returns next available B<TokenNumber> and any matched B<TokenText> from
input stream by removing it from the input stream. Token number and text of
zero corresponds to end of input (EOI).

=item B<Peek>

    ($TokenNumber, $TokenText) = $YYLexer->Peek();

Returns next available B<TokenNumber> and any matched B<TokenText> from
input stream by simply looking ahead and without removing it from the input stream.
Token number and text of zero corresponds to end of input (EOI).

=item B<SetupYYTabFile>

    $YYLexer = $YYLexer->SetupYYTabFile($YYTabFile);

Processes token labels to integers data map in specified I<YYTabFile> and returns
I<YYLexer>.

Notes:

    . YYTabFile must be a complete path or available through @INC path in the
      same directory where this package is located.
    . Name of YYTabFile might start with any valid sub directory name in @INC
      For example, "Parsers/<YYTablFile>" implies the tab file in parsers sub
      directory under MayaChemTools lib directory as it would be already in @INC
      path.
    . YYTabFile must be explicitly set by the caller. The default YYTabFile name,
      y.tab.ph, generated by byacc is not used implicitly to avoid confusion
      among multiple distinct instances of YYLexer.
    . YYTabFile is generated by byacc during its usage with -d options and
      contains mapping of token codes to token names/labels. YYLexer used this
      file to map token labels to token codes before returning token code and
      value pair back to yyparse function used by byacc.
    . User defined token numbers start from 257
    . Token number for any user defined token EOI is mapped to its value before
      default token number of 0 for EOI.

    The format of YYTabFile generated by byacc during generation of parser code in
    Perl code is:

    ... ...
    $NUMBER=257;
    $ADDOP=258;
    $SUBOP=259;
    ... ..

=item B<YYLex>

    ($TokenNumber, $TokenText) = $YYLexer->YYLex();
    ($TokenNumber, $TokenText) = $YYLexer->YYLex($Mode);

Returns available B<TokenNumber> and any matched B<TokenText> from
input stream by either removing it from the input stream or by simply looking
ahead and without removing it from the input stream. Token number and text of
zero corresponds to end of input (EOI).

Possible I<Mode> values: I<Peek, Next>. Default: I<Next>.

I<YYLex> is designed to be used with B<yyparse> code generated by running
B<byacc> on a parsers defined using parser definition B<ParserName.yy> file.

Notes:

    . Token label and value pairs returned by Lexer from by base class, which
       can't be mapped to token labels specified in YYTabFile are ignored.
    . Token text of length 1 returned by Lexer from base class without a
       corresponding explicit token label, which can't be mapped to a token
       number using Perl ord function, is ignored.

=item B<StringifyYYLexer>

    $YYLexerString = $YYLexer->StringifyYYLexer();

Returns a string containing information about I<YYLexer> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Lexer.pm, SimpleCalcYYLexer.pm, SimpleCalcParser.yy

=head1 COPYRIGHT

Copyright (C) 2017 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
