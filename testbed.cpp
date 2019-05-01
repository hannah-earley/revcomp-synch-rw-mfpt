#include <stack>

#if defined _NSIG
  #define NUM_SIGS _NSIG
#elif defined SIGRTMAX
  #define NUM_SIGS SIGRTMAX
#else
  #define NUM_SIGS 32
#endif

static std::stack<struct sigaction> oldactstack[NUM_SIGS];

void ignore_once(int signum) {
    char buf[11]="interrupt\n";
    if (write(1, buf, 10) != 10) {}
    if (signum < NUM_SIGS && !oldactstack[signum].empty()) {
        struct sigaction oldact = oldactstack[signum].top();
        sigaction(signum, &oldact, NULL);
        oldactstack[signum].pop();
    }
}

void handler(int signum) {
    if (write(1, "hello\n", 6) != 6) {}
}

int install_ign_once(int signum) {
    struct sigaction newact, oldact;
    newact.sa_handler = ignore_once;
    sigemptyset (&newact.sa_mask);
    newact.sa_flags = SA_RESETHAND;
    int ret = sigaction(signum, &newact, &oldact);
    if (0 <= signum && signum < NUM_SIGS)
        oldactstack[signum].push(oldact);
    return ret;
}

int sigpoll1(int signal) {
    sigset_t sigpend, sigmask;
    sigemptyset(&sigmask);
    sigaddset(&sigmask,signal);
    sigpending(&sigpend);
    if (sigismember(&sigpend, signal)) {
        int sigret;
        if (sigwait(&sigmask,&sigret) == 0)
            return sigret;
    }
    return -1;
}

int sigpoll(sigset_t *sigmask) {
    // sigtimedwait not available on macos
    // here we simulate its polling function...
    sigset_t sigpend;
    sigpending(&sigpend);
    for (int sig = 0; sig < NUM_SIGS; sig++) {
        if (sigismember(sigmask, sig) && sigismember(&sigpend, sig)) {
            int sigret;
            if (sigwait(sigmask,&sigret) == 0)
                return sigret;
            break;
        }
    }
    return -1;
}

void testbed() {
    signal(SIGINT, handler);
    sigset_t sigmask,oldmask; //,sigpend;
    // struct timespec poll = {0,0};

    sigemptyset(&sigmask);
    sigaddset(&sigmask,SIGINT);
    sigaddset(&sigmask,SIGTERM);
    sigprocmask(SIG_BLOCK,&sigmask,&oldmask);

    bool cont = true;
    #pragma omp parallel
    {
        volatile size_t _ = 0;
        int thread = omp_get_thread_num();
        std::string out = std::to_string(thread) + "\n";

        while (cont) {
            for (_ = 1000000000; _ > 0; _--);
            std::cout << out;

            if (thread == 0) {
                // sigpending(&sigpend);
                // if (sigismember(&sigpend, SIGINT)) {
                //     printf("received ctrl-c\n");
                //     cont = false;
                // }
                // if (sigismember(&sigpend, SIGTERM)) {
                //     printf("received sigterm\n");
                //     cont = false;
                // }
                int caught = sigpoll(&sigmask);
                if (caught >= 0) {
                    printf("caught: %d\n", caught);
                    cont = false;
                }
            }
        }
    }

    install_ign_once(SIGINT);
    install_ign_once(SIGTERM);
    sigprocmask(SIG_SETMASK,&oldmask,NULL);
    printf("returning from func\n");

    for(int i=0;i<10;i++) {
        sleep(1);
        printf("%d\n",i);
    }
    return;
}
