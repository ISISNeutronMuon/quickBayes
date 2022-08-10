### Windows and MinGW

We currently rely on an external fortran compiler, `tdm64-gcc 4.6.1`, as the current code is sensitive
to the compiler version.
This compiler has come from [tdm64-gcc-4.6.1.exe](https://sourceforge.net/projects/tdm-gcc/files/TDM-GCC%20Installer/Previous/1.1006.0/tdm64-gcc-4.6.1.exe/download) and [gcc-4.6.1-tdm64-1-fortran.zip](https://sourceforge.net/projects/tdm-gcc/files/TDM-GCC%20Old%20Releases/TDM-GCC%204.6%20series/4.6.1-tdm64-1/gcc-4.6.1-tdm64-1-fortran.zip/download) and has been packaged together in a single zxip file for convenience.
To install:

- Download MinGW64.7z and unzip it into ``C:\MinGW64``
- Add ``C:\MinGW64\bin`` to your ``PATH`` environment variable ([instructions here](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/))
- Restart any terminal or powershell instances to capture the new environment variable settings

