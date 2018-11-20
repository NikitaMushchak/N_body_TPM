#define macros
EXECUTABLE_NAME = tpm.exe
CPP_FLAGS = /Ox /favor:%PROCESSOR_ARCHITECTURE% \
	/D_TARGET_="%PROCESSOR_ARCHITECTURE%" \
	/D_TIMESTAMP_="%date% %time%" \
	/D_PLATFORM_="%OS%"
DIR_SRC = .\source
DIR_DATA = $(DIR_SRC)\data
DIR_BUILD = .\build
SRC_FILES= $(DIR_SRC)\*.cpp
DATA_FILES= $(DIR_DATA)\*.dat

#SHELL = powershell.exe -NoLogo -NoProfile -NonInteractive -ExecutionPolicy Bypass -Command

all: clean create_dirs tpm



clean:
	@if exist $(DIR_BUILD) rmdir /S /Q $(DIR_BUILD)

create_dirs:
	@if not exist $(DIR_BUILD) mkdir $(DIR_BUILD)
	@if not exist $(DIR_BUILD)\InitialConditions mkdir $(DIR_BUILD)\InitialConditions
	@if not exist $(DIR_BUILD)\Results mkdir $(DIR_BUILD)\Results
	@if not exist $(DIR_BUILD)\Results\Concentration mkdir $(DIR_BUILD)\Results\Concentration
	@if not exist $(DIR_BUILD)\Results\Opening mkdir $(DIR_BUILD)\Results\Opening
	@if exist $(DIR_DATA) if exist $(DATA_FILES) copy $(DATA_FILES) $(DIR_BUILD)\InitialConditions

$(EXECUTABLE_NAME) : $(SRC_FILES)
	cl /EHsc $(CPP_FLAGS) /Fe$(DIR_BUILD)\$(EXECUTABLE_NAME) /I$(DIR_INCLUDE) $(SRC_FILES)
	del *.obj


tpm: $(EXECUTABLE_NAME)
