package FileIO::FileIO;
#
# $RCSfile: FileIO.pm,v $
# $Date: 2017/01/12 18:59:14 $
# $Revision: 1.29 $
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
use FileHandle;
use ObjectProperty;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(ObjectProperty Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeFileIO();

  $This->_InitializeFileIOProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeFileIO {
  my($This) = @_;

  # File name...
  $This->{Name} = '';

  # Read, write or append...
  $This->{Mode} = 'Read';

  # Open/close status...
  $This->{Status} = 0;

  # File handle returned by file open...
  $This->{FileHandle} = '';
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

}

# Initialize object properties....
sub _InitializeFileIOProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  return $This;
}

# Close any open file...
sub DESTROY {
  my($This) = @_;

  $This->Close();

  return $This;
}

# Set file name and make sure it's not already set...
#
sub SetName {
  my($This, $Name) = @_;

  if ($This->{Name}) {
    croak "Error: ${ClassName}->SetName: Can't set file name to $Name:  $This->{Name}...";
  }

  $This->{Name} = $Name;

  return $This;
}

# Open file using specified mode...
#
sub Open {
  my($This, $Mode) = @_;

  if ($This->{Status}) {
    croak "Error: ${ClassName}->Open: Can't open file $This->{Name}: It's already open...";
  }

  if (defined $Mode) {
    # Set mode...
    $This->SetMode($Mode);
  }

  # Get name and mode...
  my($Name);
  $Name = $This->{Name};
  $Mode = $This->_GetOpenMode();

  # Open the file using specified mode and store FileHandle...
  my($FileHandle);
  $FileHandle = new FileHandle("${Mode}${Name}");
  if (!defined $FileHandle) {
    croak "Error: ${ClassName}->Open: Can't open $Name:  $! ...";
  }
  $This->{FileHandle} = $FileHandle;
  $This->{Status} = 1;

  return $This;
}

# Close an open file...
sub Close {
  my($This) = @_;

  if ($This->{Status}) {
    $This->{FileHandle}->close();
  }
  $This->{Status} = 0;

  return $This;
}

# Supported Mode values are: Read, Write, Append, <, >, >>, r, w, a
#
sub SetMode {
  my($This, $SpecifiedMode) = @_;
  my($Mode);

  if (!defined $SpecifiedMode) {
    $SpecifiedMode = 'Read';
  }

  MODE: {
    if ($SpecifiedMode =~ /^(Read|<|r)$/i) { $Mode = 'Read'; last MODE; }
    if ($SpecifiedMode =~ /^(Write|>|w)$/i) { $Mode = 'Write'; last MODE; }
    if ($SpecifiedMode =~ /^(Append|>>|a)$/i) { $Mode = 'Append'; last MODE; }
    croak "Error: ${ClassName}->SetMode: Specified mode value, $SpecifiedMode, is not valid: Supported values: Read, Write, Append, <, >, >>, r, w, a...";
  }
  $This->{Mode} = $Mode;

  return $This;
}

# Get mode values to be used for file open function: <, >, >>
#
sub _GetOpenMode {
  my($This) = @_;
  my($Mode);

  MODE: {
    if ($This->{Mode} =~ /^(Read|<|r)$/i) { $Mode = '<'; last MODE; }
    if ($This->{Mode} =~ /^(Write|>|w)$/i) { $Mode = '>'; last MODE; }
    if ($This->{Mode} =~ /^(Append|>>|a)$/i) { $Mode = '>>'; last MODE; }
    $Mode = '';
  }
  return $Mode;
}

1;

__END__

=head1 NAME

FileIO

=head1 SYNOPSIS

use FileIO::FileIO;

use FileIO::FileIO qw(:all);

=head1 DESCRIPTION

B<FIleIO> class provides following methods:

new, Close, Open, SetMode

B<FleIO> class serves as a base class for all classes involved in file IO. It is derived from
B<ObjectProperty> base class which provides methods not explicitly defined in B<Atom> or
B<ObjectProperty> class using Perl's AUTOLOAD functionality. These methods
are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

=head2 METHODS

=over 4

=item B<new>

    $NewFileIO = new FileIO(%NamesAndValues);

Using specified I<FileIO> property names and values hash, B<new> method creates a new object
and returns a reference to a newly created B<FileIO> object. By default, the following properties are
initialized:

    Name = '';
    Mode = 'Read';
    Status = 0;
    FileHandle = '';

=item B<Close>

    $FileIO->Close();

Close open file and returns I<FileIO>.

=item B<Open>

    $FileIO->Open();

Opens the file using file I<Name> and I<Mode> and returns I<FileIO>.

=item B<SetMode>

    $FileIO->SetMode($Mode);

Sets up file I<Mode> and returns I<FileIO> Default I<Mode> value: I<Read>.
Supported I<Mode> values:

    Read, Write, Append, <, >, >>, r, w, a

=item B<SetName>

    $FileIO->SetName($Name);

Sets up file I<Name> and returns I<FileIO>.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MoleculeFileIO.pm, MDLMolFileIO.pm, SDFileIO.pm

=head1 COPYRIGHT

Copyright (C) 2017 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
