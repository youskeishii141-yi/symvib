; Inno Setup script for the SymVib Windows installer (SymVib-Setup.exe).
;
; The installer copies the PyInstaller folder (dist\SymVib) itself, so users do not have to
; extract a zip with thousands of files. By default it installs for the current user only
; (no administrator rights needed, into %LOCALAPPDATA%\Programs\SymVib); the user can choose
; to install for all users instead.
;
; Build:  ISCC.exe /DAppVersion=2.2.3 installer\SymVib.iss    (output: dist\SymVib-Setup.exe)

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#ifndef SourceDir
  #define SourceDir "..\dist\SymVib"
#endif

[Setup]
AppId={{6445B75E-8E7F-4DB3-8E97-7BC4C2B30531}
AppName=SymVib
AppVersion={#AppVersion}
AppVerName=SymVib {#AppVersion}
AppPublisher=Yosuke Ishii
AppPublisherURL=https://github.com/youskeishii141-yi/symvib
AppSupportURL=https://github.com/youskeishii141-yi/symvib/issues
AppUpdatesURL=https://github.com/youskeishii141-yi/symvib/releases
DefaultDirName={autopf}\SymVib
DefaultGroupName=SymVib
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
LicenseFile={#SourceDir}\LICENSE
OutputDir=..\dist
OutputBaseFilename=SymVib-Setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
UninstallDisplayIcon={app}\SymVib.exe
UninstallDisplayName=SymVib {#AppVersion}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "japanese"; MessagesFile: "compiler:Languages\Japanese.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[InstallDelete]
; Remove the files of a previous version before copying the new ones.
Type: filesandordirs; Name: "{app}\_internal"

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\SymVib"; Filename: "{app}\SymVib.exe"
Name: "{autodesktop}\SymVib"; Filename: "{app}\SymVib.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\SymVib.exe"; Description: "{cm:LaunchProgram,SymVib}"; Flags: nowait postinstall skipifsilent
