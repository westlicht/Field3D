//-*****************************************************************************
//
// Copyright (c) 2013,
//  Sony Pictures Imageworks Inc. and
//  Industrial Light & Magic, a division of Lucasfilm Entertainment Company Ltd.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// *       Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// *       Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
// *       Neither the name of Industrial Light & Magic nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//-*****************************************************************************

#include "OStream.h"
#include "Encryption.h"

#include <themispp/secure_cell.hpp>

#include <fstream>
#include <stdexcept>

namespace Alembic {
namespace Ogawa {
namespace ALEMBIC_VERSION_NS {

class OStream::PrivateData
{
public:
    PrivateData(const std::string & iFileName) :
        stream(NULL), fileName(iFileName), startPos(0)
    {
        std::ofstream * filestream = new std::ofstream(fileName.c_str(),
            std::ios_base::trunc | std::ios_base::binary);
        if (filestream->is_open())
        {
            stream = filestream;
            stream->exceptions ( std::ofstream::failbit |
                                 std::ofstream::badbit );
        }
        else
        {
            filestream->close();
            delete filestream;
        }
    }

    PrivateData(std::ostream * iStream) : stream(iStream), startPos(0)
    {
        if (stream)
        {
            stream->exceptions ( std::ostream::failbit |
                                 std::ostream::badbit );

            startPos = stream->tellp();
            if (startPos == INVALID_DATA)
            {
                throw std::runtime_error("Illegal start of Ogawa stream");
            }
        }
    }

    ~PrivateData()
    {
        // if this was done via file, try to clean it up
        if (!fileName.empty() && stream)
        {
            std::ofstream * filestream = dynamic_cast<std::ofstream *>(stream);
            if (filestream)
            {
                filestream->close();
                delete filestream;
            }
        }
    }

    std::ostream * stream;
    std::string fileName;
    Alembic::Util::uint64_t startPos;
    Alembic::Util::mutex lock;
};

OStream::OStream(const std::string & iFileName) :
    mData(new PrivateData(iFileName))
{
    init();
}

// we'll be writing from this already open stream which we don't own
OStream::OStream(std::ostream * iStream) : mData(new PrivateData(iStream))
{
    init();
}

OStream::~OStream()
{
    // write our "frozen" byte (totally done writing)
    if (isValid())
    {
        char frozen = 0xff;
        mData->stream->seekp(mData->startPos + 5).write(&frozen, 1).flush();
    }

    // Encrypt file if requested
    if (isValid() && Encryption::encryptSettings.enabled) {
        // close plain file
        std::string fileName = mData->fileName;
        mData.reset();

        // read plain file to memory
        std::ifstream ifs(fileName.c_str(), std::ios::in | std::ios::binary);
        auto p = ifs.tellg();
        ifs.seekg(0, std::ios::end);
        size_t size = ifs.tellg() - p;
        ifs.seekg(p);
        std::vector<uint8_t> plain(size);
        ifs.read(reinterpret_cast<char *>(plain.data()), plain.size());
        ifs.close();

        // convert password
        std::vector<uint8_t> passwordData;
        const auto &password = Encryption::encryptSettings.password;
        std::copy(password.begin(), password.end(), std::back_inserter(passwordData));

        // encrypt data
        themispp::secure_cell_seal_t sm(passwordData);
        auto encrypted = sm.encrypt(plain);

        // write encrypted file
        std::ofstream ofs(fileName.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        ofs.write(reinterpret_cast<const char *>(encrypted.data()), encrypted.size());
        ofs.close();
    }    
}

bool OStream::isValid()
{
    return mData->stream != NULL;
}

void OStream::init()
{
    // simple temporary endian check
    union {
        Util::uint32_t l;
        char c[4];
    } u;

    u.l = 0x01234567;

    if (u.c[0] != 0x67)
    {
        throw std::runtime_error(
            "Ogawa currently only supports little-endian writing.");
    }

    if (isValid())
    {
        const char header[] = {
            'O', 'g', 'a', 'w', 'a',  // special magic number
            0,       // this will be 0xff when the entire archive is done
            0, 1,    // 16 bit format version number
            0, 0, 0, 0, 0, 0, 0, 0}; // position of the first group
        mData->stream->write(header, sizeof(header)).flush();
    }
}

Alembic::Util::uint64_t OStream::getAndSeekEndPos()
{
    if (isValid())
    {
        Alembic::Util::scoped_lock l(mData->lock);
        Alembic::Util::uint64_t lastp =
            mData->stream->seekp(0, std::ios_base::end).tellp();
        if (lastp == INVALID_DATA || lastp < mData->startPos)
        {
            throw std::runtime_error(
                "Illegal position returned Ogawa::OStream::getAndSeekEndPos");

            return 0;
        }
        return lastp - mData->startPos;
    }
    return 0;
}

void OStream::seek(Alembic::Util::uint64_t iPos)
{
    if (isValid())
    {
        Alembic::Util::scoped_lock l(mData->lock);
        mData->stream->seekp(iPos + mData->startPos);
    }
}

void OStream::write(const void * iBuf, Alembic::Util::uint64_t iSize)
{
    if (isValid())
    {
        Alembic::Util::scoped_lock l(mData->lock);
        mData->stream->write((const char *)iBuf, iSize).flush();
    }
}

} // End namespace ALEMBIC_VERSION_NS
} // End namespace Ogawa
} // End namespace Alembic
