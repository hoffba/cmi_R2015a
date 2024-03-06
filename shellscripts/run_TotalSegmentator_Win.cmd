:: Inputs:
::	1: Image filename
::	2: Saved segmentation filename

:: First make sure python is set up to run TotalSegmentator
SET stat=false
for /f %%i in ('conda env list') do (
	if "%%i" == "TSenv" if not %stat% == true (SET stat=true & goto ACTIVATE_TSenv)
)

:ACTIVATE_TSenv
echo Setting up environment
SET ACTIVATEPATH=C:\ProgramData\anaconda3\Scripts\activate.bat
SET TSpath=%userprofile%\.conda\envs\TSenv
if %stat% == true (
	call %ACTIVATEPATH% TSenv
) else (
	conda create -y -n TSenv
	call %ACTIVATEPATH% TSenv
	conda install -y pytorch torchvision
	pip install totalsegmentator
)

:: Find TotalSegmentator.py full path
echo Finding TotalSegmentator
for /f %%i in ('dir /a:-d /b/s %userprofile%\.conda\envs\TSenv\*TotalSegmentator.py') do (
	SET TSpath=%%i
	goto RUN_TS
)

:RUN_TS
echo Running TotalSegmentator
:: python %TSpath% -i %1 -o %2 --ml
echo %TSpath%