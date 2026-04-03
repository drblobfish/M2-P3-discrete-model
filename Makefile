main:
	g++ main.cpp CellularPotts.cpp -o main.exe -ggdb

UI:
	g++ -Iextern/ -Iextern/raylib/include/ UI.cpp CellularPotts.cpp extern/raylib/lib/libraylib.a -lGL -lm -lpthread -ldl -lrt -lX11 -o UI.exe -ggdb
