#ifndef _Alembic_Ogawa_Encryption_h_
#define _Alembic_Ogawa_Encryption_h_

#include <string>

class Encryption {
public:
    struct Settings {
        bool enabled = false;
        std::string password;
    };

    static Settings encryptSettings;
    static Settings decryptSettings;
};


#endif // _Alembic_Ogawa_Encryption_h_
