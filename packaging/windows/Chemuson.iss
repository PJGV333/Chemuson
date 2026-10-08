; Instalador MVP de Chemuson para Inno Setup.

#define MyAppName "ChemUSON"
#define MyAppVersion GetEnv("CHEMUSON_VERSION")
#define MyAppPublisher "ChemUSON"
#define MyAppExeName "Chemuson.exe"
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
OutputDir=..\..\dist-installer
OutputBaseFilename=Chemuson-v{#MyAppVersion}-windows-x86_64-setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "spanish"; MessagesFile: "compiler:Languages\\Spanish.isl"

[Files]
Source: "..\\..\\dist\\{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
Source: "chemuson-installed.marker"; DestDir: "{app}"; DestName: ".chemuson-installed"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\\ChemUSON"; Filename: "{app}\\{#MyAppExeName}"
Name: "{autodesktop}\\ChemUSON"; Filename: "{app}\\{#MyAppExeName}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Crear acceso directo en escritorio"; GroupDescription: "Opciones adicionales:"

[Run]
Filename: "{app}\\{#MyAppExeName}"; Description: "Iniciar ChemUSON"; Flags: nowait postinstall skipifsilent
