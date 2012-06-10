#
# Note: popt is now an rpm sub-package (including libpopt.so*) so you probably
# shouldn't need to use this spec file to package popt anymore.
#
Summary: A C library for parsing command line parameters.
Name: popt
Version: 1.10.2
Release: 0.1
License: X Consortium
Group: System Environment/Libraries
Source: ftp://ftp.redhat.com/pub/redhat/code/popt/popt-%{version}.tar.gz
BuildRoot: /var/tmp/%{name}root

%description
Popt is a C library for parsing command line parameters.  Popt
was heavily influenced by the getopt() and getopt_long() functions,
but it improves on them by allowing more powerful argument expansion.
Popt can parse arbitrary argv[] style arrays and automatically set
variables based on command line arguments.  Popt allows command
line arguments to be aliased via configuration files and includes
utility functions for parsing arbitrary strings into argv[] arrays
using shell-like rules.

Install popt if you're a C programmer and you'd like to use its
capabilities.

%prep
%setup -q

%build
#CFLAGS="$RPM_OPT_FLAGS" ./configure --prefix=/usr

%configure
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{_prefix}/lib/libpopt.*
%{_prefix}/include/popt.h
%{_prefix}/man/man3/popt.3
%{_prefix}/share/locale/*/LC_MESSAGES/popt.mo

%changelog
* Thu Dec 10 1998 Michael Johnson <johnsonm@redhat.com>
- released 1.2.2; see CHANGES

* Tue Nov 17 1998 Michael K. Johnson <johnsonm@redhat.com>
- added man page to default install

* Thu Oct 22 1998 Erik Troan <ewt@redhat.com>
- see CHANGES file for 1.2

* Thu Apr 09 1998 Erik Troan <ewt@redhat.com>
- added ./configure step to spec file
