/*
Copyright (c) 2018 Xiao Wang (wangxiao@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Enquiries about further applications and development opportunities are welcome.
*/

#include "utils/io.h"

NetIO::NetIO(const char* address, int port, bool full_buffer, bool quiet) {
    this->port = port;
    is_server  = (address == nullptr);
    if (address == nullptr) {
        struct sockaddr_in dest;
        struct sockaddr_in serv;
        socklen_t socksize = sizeof(struct sockaddr_in);
        memset(&serv, 0, sizeof(serv));
        serv.sin_family      = AF_INET;
        serv.sin_addr.s_addr = htonl(INADDR_ANY); /* set our address to any interface */
        serv.sin_port        = htons(port);       /* set the server port number */
        mysocket             = socket(AF_INET, SOCK_STREAM, 0);
        int reuse            = 1;
        setsockopt(mysocket, SOL_SOCKET, SO_REUSEADDR, (const char*)&reuse, sizeof(reuse));
        if (::bind(mysocket, (struct sockaddr*)&serv, sizeof(struct sockaddr)) < 0) {
            perror("error: bind");
            exit(1);
        }
        if (listen(mysocket, 1) < 0) {
            perror("error: listen");
            exit(1);
        }
        consocket = accept(mysocket, (struct sockaddr*)&dest, &socksize);
        close(mysocket);
    }
    else {
        addr = string(address);

        struct sockaddr_in dest;
        memset(&dest, 0, sizeof(dest));
        dest.sin_family      = AF_INET;
        dest.sin_addr.s_addr = inet_addr(address);
        dest.sin_port        = htons(port);

        while (1) {
            consocket = socket(AF_INET, SOCK_STREAM, 0);

            if (connect(consocket, (struct sockaddr*)&dest, sizeof(struct sockaddr)) == 0) {
                break;
            }

            close(consocket);
            usleep(1000);
        }
    }
    set_nodelay();
    stream = fdopen(consocket, "wb+");
    buffer = new char[NETWORK_BUFFER_SIZE];
    memset(buffer, 0, NETWORK_BUFFER_SIZE);
    if (full_buffer) {
        setvbuf(stream, buffer, _IOFBF, NETWORK_BUFFER_SIZE);
    }
    else {
        setvbuf(stream, buffer, _IONBF, NETWORK_BUFFER_SIZE);
    }
    this->FBF_mode = full_buffer;
    if (!quiet)
        cout << "connected\n";
}

NetIO::~NetIO() {
    fflush(stream);
    close(consocket);
    delete[] buffer;
}

void NetIO::sync() {
    int tmp = 0;
    if (is_server) {
        send_data(&tmp, 1);
        recv_data(&tmp, 1);
    }
    else {
        recv_data(&tmp, 1);
        send_data(&tmp, 1);
        flush();
    }
}

void NetIO::send_data(const void* data, int len, bool count_comm) {
    if (count_comm) {
        if (last_call != LastCall::Send) {
            num_rounds++;
            last_call = LastCall::Send;
        }
        counter += len;
    }
    int sent = 0;
    while (sent < len) {
        int res = fwrite(sent + (char*)data, 1, len - sent, stream);
        if (res >= 0) {
            sent += res;
        }
        else {
            fprintf(stderr, "error: net_send_data %d\n", res);
        }
    }
    has_sent = true;
}

void NetIO::recv_data(void* data, int len, bool count_comm) {
    if (count_comm) {
        if (last_call != LastCall::Recv) {
            num_rounds++;
            last_call = LastCall::Recv;
        }
    }
    if (has_sent) {
        fflush(stream);
    }
    has_sent = false;
    int sent = 0;
    while (sent < len) {
        int res = fread(sent + (char*)data, 1, len - sent, stream);
        if (res >= 0) {
            sent += res;
        }
        else {
            fprintf(stderr, "error: net_send_data %d\n", res);
        }
    }
}

IOPack::IOPack(int party, string address) {
    if (party == ALICE) {
        this->io     = new NetIO(nullptr, ALICE_SEND_PORT, false, true);
        this->io_rev = new NetIO(nullptr, BOB_SEND_PORT);
    }
    else {
        this->io_rev = new NetIO(address.c_str(), ALICE_SEND_PORT, false, true);
        this->io     = new NetIO(address.c_str(), BOB_SEND_PORT);
    }
}

IOPack::~IOPack() {
    delete io;
    delete io_rev;
}

void IOPack::send_data(const void* data, int len, bool count_comm) const {
    io->send_data(data, len);
    io->last_call = LastCall::Send;
    io->last_call = LastCall::Send;
}

void IOPack::recv_data(void* data, int len, bool count_comm) const {
    io_rev->recv_data(data, len);
    io->last_call = LastCall::Recv;
    io->last_call = LastCall::Recv;
}