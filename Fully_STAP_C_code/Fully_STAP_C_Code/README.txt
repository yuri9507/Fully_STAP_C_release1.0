Project : Fully_STAP_C_Code
pre-install : 1. GSL 2.Visual Studio 2019 3.git
------------------------------------------------
1.Installing GSL on Windows / Linux 
@Source : https://solarianprogrammer.com/2020/01/26/getting-started-gsl-gnu-scientific-library-windows-macos-linux/
Open a Windows "PowerShell window"
cd C:\
mkdir DEV
cd DEV
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg integrate install
.\vcpkg install gsl gsl:x64-windows
 
2.open the Fully_STAP_C_Code.sln with Visual Studio 2019 
3.run it 