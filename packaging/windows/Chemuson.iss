; Instalador MVP de Chemuson para Inno Setup.

#define MyAppName "ChemUSON"
#define MyAppVersion GetEnv("CHEMUSON_VERSION")
#define MyAppPublisher "ChemUSON"
#define MyAppExeName "Chemuson.exe"
#define BuildOutputRoot GetEnv("CHEMUSON_OUTPUT_DIR")
#if MyAppVersion == ""
  #error "CHEMUSON_VERSION must be set explicitly before compiling the installer."
#endif

[Setup]
AppId={{E2D4477C-AE35-4C8E-9F7E-8C2E4DBE69A7}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
UninstallDisplayIcon={app}\{#MyAppExeName}
AppMutex=ChemusonMainMutex
DefaultDirName={autopf}\Chemuson
DefaultGroupName=ChemUSON
DisableProgramGroupPage=yes
#if BuildOutputRoot == ""
OutputDir=..\..\dist-installer
#else
OutputDir={#BuildOutputRoot}\dist-installer
#endif
OutputBaseFilename=Chemuson-v{#MyAppVersion}-windows-x86_64-setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern
SetupIconFile=..\..\assets\branding\chemuson.ico

[Languages]
Name: "spanish"; MessagesFile: "compiler:Languages\\Spanish.isl"

[Files]
#if BuildOutputRoot == ""
Source: "..\\..\\dist\\{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
#else
Source: "{#BuildOutputRoot}\\dist\\{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
#endif
Source: "chemuson-installed.marker"; DestDir: "{app}"; DestName: ".chemuson-installed"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\\ChemUSON"; Filename: "{app}\\{#MyAppExeName}"
Name: "{autodesktop}\\ChemUSON"; Filename: "{app}\\{#MyAppExeName}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Crear acceso directo en escritorio"; GroupDescription: "Opciones adicionales:"

[Run]
Filename: "{app}\\{#MyAppExeName}"; Description: "Iniciar ChemUSON"; Flags: nowait postinstall skipifsilent
