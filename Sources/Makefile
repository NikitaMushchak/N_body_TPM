EXECUTABLE_NAME = tpm.exe
CPP_FLAGS = /Ox /favor:%PROCESSOR_ARCHITECTURE% \
	/D_TARGET_="%PROCESSOR_ARCHITECTURE%" \
	/D_TIMESTAMP_="%date% %time%" \
	/D_PLATFORM_="win"
DIR_SRC = .\source
DIR_DATA = $(DIR_SRC)\data
DIR_BUILD = .\build
SRC_FILES = $(DIR_SRC)\*.cpp


#SHELL = powershell.exe -NoLogo -NoProfile -NonInteractive -ExecutionPolicy Bypass -Command

all: clean create_dirs tpm



clean:
	@if exist $(DIR_BUILD) rmdir /S /Q $(DIR_BUILD)

create_dirs:
	@if not exist $(DIR_BUILD) mkdir $(DIR_BUILD)
	@if not exist $(DIR_BUILD)\Results mkdir $(DIR_BUILD)\Results




$(EXECUTABLE_NAME) : $(SRC_FILES)
	cl /EHsc $(CPP_FLAGS) /Fe$(DIR_BUILD)\$(EXECUTABLE_NAME) $(SRC_FILES)
	del *.obj


tpm: $(EXECUTABLE_NAME)
