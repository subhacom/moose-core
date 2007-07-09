#include <iostream>
#include <string>
#include <vector>
#include<stdio.h>
#include<termios.h>
#include<unistd.h>
#include<sys/select.h>
#include<sys/time.h>

using namespace std;

//kbhit, Non-blocking keypress detector, when go keypress, return 1 else always return 0
int kbhit()
{
    struct timeval tv;
    fd_set fds;
    tv.tv_sec = 0;
    tv.tv_usec = 0;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds); //STDOUT_FILENO is 0
    select(STDIN_FILENO+1, &fds, NULL, NULL, &tv);
    return FD_ISSET(STDIN_FILENO, &fds);
}

#define NB_DISABLE 0
#define NB_ENABLE 1

// This function changes the terminal state. It turns out that
// the canonical state is highly desirable. So I don't actually
// use this function at this time.
// The only missing link here is that I still cannot trap control-p.
void nonblock(int state)
{
    struct termios ttystate;

    //get the terminal state
    tcgetattr(STDIN_FILENO, &ttystate);

    if (state==NB_ENABLE)
    {
        //turn off canonical mode
        ttystate.c_lflag &= ~ICANON;
        //minimum of number input read.
        ttystate.c_cc[VMIN] = 1; 
    }
    else if (state==NB_DISABLE)
    {
        //turn on canonical mode
        ttystate.c_lflag |= ICANON;
    }
    //set the terminal attributes.
    tcsetattr(STDIN_FILENO, TCSANOW, &ttystate);

}

bool nonBlockingGetLine( string& s )
{
	static vector< string > history;
	static unsigned int historyIndex = 0;
	static char line[400];

	usleep( 10 );
	if ( kbhit() ) {
		fgets( line, 399, stdin );
		s = line;
		history.push_back( s );
		return 1;

		/*
		cout << "." << flush;
		char c = fgetc( stdin );
		switch ( c ) {
			case 0x10: // control-p
				break;
			case 0x0e: // control-n
				break;
			case '\b': // backspace
			case 0x7f: // delete
				// do something!
//				fputc( c, stdout );
				break;
			case 0x15: // clear line, hopefully control-u
				break;
			case '\n':
			case '\r':
				if ( s.find_first_not_of( " \t" ) != s.npos ) {
					history.push_back( s );
				}
				s.push_back( c );
//				fputc( c, stdout );
				return 1;
				break;
			default:
//				fputc( c, stdout );
				s.push_back( c );
				break;
		}
		*/
	}
	return 0;
}
