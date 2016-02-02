OMBlast
================
An alignment tool for optical mapping data

Please refer to "OMBlastManual.pdf" for more details. 

Quick steps
------------
1.	Compile the OMBlast package in the OMBlast folder:
javac -d bin -sourcepath src -cp "lib/*" @classes
2.	Build a runnable jar file for OMBlast:
jar cvfm OMBlast.jar manifest -C bin .
3.	Run OMBlast:
java -jar OMBlast.jar
4.	Run other scripts:
(Linux) 		java -cp “bin:lib/*” aldenjava.script.ScriptName
(Windows) 	java -cp “bin;lib/*” aldenjava.script.ScriptName

