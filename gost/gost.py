'''
Created on Jan 7, 2014

@author: Stefan Koelbl
'''

import random

class GOST:
    """
    A basic implementation of the round function of the
    GOST hash function.
    """
    # Constant values for KeySchedule function
    C = [
        0xb1085bda1ecadae9ebcb2f81c0657c1f2f6a76432e45d016714eb88d7585c4fc4b7ce09192676901a2422a08a460d31505767436cc744d23dd806559f2a64507,
        0x6fa3b58aa99d2f1a4fe39d460f70b5d7f3feea720a232b9861d55e0f16b501319ab5176b12d699585cb561c2db0aa7ca55dda21bd7cbcd56e679047021b19bb7,
        0xf574dcac2bce2fc70a39fc286a3d843506f15e5f529c1f8bf2ea7514b1297b7bd3e20fe490359eb1c1c93a376062db09c2b6f443867adb31991e96f50aba0ab2,
        0xef1fdfb3e81566d2f948e1a05d71e4dd488e857e335c3c7d9d721cad685e353fa9d72c82ed03d675d8b71333935203be3453eaa193e837f1220cbebc84e3d12e,
        0x4bea6bacad4747999a3f410c6ca923637f151c1f1686104a359e35d7800fffbdbfcd1747253af5a3dfff00b723271a167a56a27ea9ea63f5601758fd7c6cfe57,
        0xae4faeae1d3ad3d96fa4c33b7a3039c02d66c4f95142a46c187f9ab49af08ec6cffaa6b71c9ab7b40af21f66c2bec6b6bf71c57236904f35fa68407a46647d6e,
        0xf4c70e16eeaac5ec51ac86febf240954399ec6c7e6bf87c9d3473e33197a93c90992abc52d822c3706476983284a05043517454ca23c4af38886564d3a14d493,
        0x9b1f5b424d93c9a703e7aa020c6e41414eb7f8719c36de1e89b4443b4ddbc49af4892bcb929b069069d18d2bd1a5c42f36acc2355951a8d9a47f0dd4bf02e71e,
        0x378f5a541631229b944c9ad8ec165fde3a7d3a1b258942243cd955b7e00d0984800a440bdbb2ceb17b2b8a9aa6079c540e38dc92cb1f2a607261445183235adb,
        0xabbedea680056f52382ae548b2e4f3f38941e71cff8a78db1fffe18a1b3361039fe76702af69334b7a1e6c303b7652f43698fad1153bb6c374b4c7fb98459ced,
        0x7bcd9ed0efc889fb3002c6cd635afe94d8fa6bbbebab076120018021148466798a1d71efea48b9caefbacd1d7d476e98dea2594ac06fd85d6bcaa4cd81f32d1b,
        0x378ee767f11631bad21380b00449b17acda43c32bcdf1d77f82012d430219f9b5d80ef9d1891cc86e71da4aa88e12852faf417d5d9b21b9948bc924af11bd720
        ]

    A = [
        0x641c314b2b8ee083, 0xc83862965601dd1b,
        0x8d70c431ac02a736, 0x07e095624504536c,
        0x0edd37c48a08a6d8, 0x1ca76e95091051ad,
        0x3853dc371220a247, 0x70a6a56e2440598e,
        0xa48b474f9ef5dc18, 0x550b8e9e21f7a530,
        0xaa16012142f35760, 0x492c024284fbaec0,
        0x9258048415eb419d, 0x39b008152acb8227,
        0x727d102a548b194e, 0xe4fa2054a80b329c,
        0xf97d86d98a327728, 0xeffa11af0964ee50,
        0xc3e9224312c8c1a0, 0x9bcf4486248d9f5d,
        0x2b838811480723ba, 0x561b0d22900e4669,
        0xac361a443d1c8cd2, 0x456c34887a3805b9,
        0x5b068c651810a89e, 0xb60c05ca30204d21,
        0x71180a8960409a42, 0xe230140fc0802984,
        0xd960281e9d1d5215, 0xafc0503c273aa42a,
        0x439da0784e745554, 0x86275df09ce8aaa8,
        0x0321658cba93c138, 0x0642ca05693b9f70,
        0x0c84890ad27623e0, 0x18150f14b9ec46dd,
        0x302a1e286fc58ca7, 0x60543c50de970553,
        0xc0a878a0a1330aa6, 0x9d4df05d5f661451,
        0xaccc9ca9328a8950, 0x4585254f64090fa0,
        0x8a174a9ec8121e5d, 0x092e94218d243cba,
        0x125c354207487869, 0x24b86a840e90f0d2,
        0x486dd4151c3dfdb9, 0x90dab52a387ae76f,
        0x46b60f011a83988e, 0x8c711e02341b2d01,
        0x05e23c0468365a02, 0x0ad97808d06cb404,
        0x14aff010bdd87508, 0x2843fd2067adea10,
        0x5086e740ce47c920, 0xa011d380818e8f40,
        0x83478b07b2468764, 0x1b8e0b0e798c13c8,
        0x3601161cf205268d, 0x6c022c38f90a4c07,
        0xd8045870ef14980e, 0xad08b0e0c3282d1c,
        0x47107ddd9b505a38, 0x8e20faa72ba0b470
        ]

    Ainverse = [
        0xdc6e36e95e0a310b,
        0xa5dc6ccfbc146216,
        0x57a5d8836528c42c,
        0xae57ad1bca509558,
        0x41ae473689a037b0,
        0x82418e6c0f5d6e7d,
        0x198201d81ebadcfa,
        0x321902ad3c69a5e9,
        0x49d9d1cb12a5f7b1,
        0x92afbf8b2457f37f,
        0x3943630b48aefbfe,
        0x7286c6169041ebe1,
        0xe411912c3d82cbdf,
        0xd5223f587a198ba3,
        0xb7447eb0f4320b5b,
        0x7388fc7df56416b6,
        0x32786fdea3a318a8,
        0x64f0dea15b5b304d,
        0xc8fda15fb6b6609a,
        0x8de75fbe7171c029,
        0x07d3be61e2e29d52,
        0x0ebb61c2d9d927a4,
        0x1c6bc299afaf4e55,
        0x38d6992f43439caa,
        0x70256d11ef564e51,
        0xe04ada22c3ac9ca2,
        0xdd94a9449b452559,
        0xa7354f882b8a4ab2,
        0x536a9e0d56099479,
        0xa6d4211aac1235f2,
        0x51b5423445246af9,
        0xa27784688a48d4ef,
        0x9813d5ff7fe3f8ad,
        0x2d26b7e3fedbed47,
        0x5a4c73dbe1abc78e,
        0xb498e6abdf4b9301,
        0x752dd14ba3963b02,
        0xea5abf965b317604,
        0xc9b46331b662ec08,
        0x8f75c66271c4c510,
        0x91073f14f6277fd0,
        0x3f0e7e28f14efebd,
        0x7e1cfc50ff9ce167,
        0xfc38e5a0e325dfce,
        0xe570d75ddb4aa381,
        0xd7e0b3baab945b1f,
        0xb3dd7b694b35b63e,
        0x7ba7f6d2966a717c,
        0xabf2017c29f19997,
        0x4bf902f852ff2f33,
        0x96ef04eda4e35e66,
        0x31c308c755dbbccc,
        0x629b1093aaab6585,
        0xc42b203b494bca17,
        0x955640769296892e,
        0x37ac80ec39310f5c,
        0xffb3d920c3520bec,
        0xe37baf409ba416c5,
        0xdbf643802b552c97,
        0xabf1861d56aa5833,
        0x4bff113aac49b066,
        0x96e3227445927dcc,
        0x31db44e88a39fa85,
        0x62ab88cd0972e917
        ]

    sub = [
            252, 238, 221, 17, 207, 110, 49, 22, 251, 196, 250,
            218, 35, 197, 4, 77, 233, 119, 240, 219, 147, 46,
            153, 186, 23, 54, 241, 187, 20, 205, 95, 193, 249,
            24, 101, 90, 226, 92, 239, 33, 129, 28, 60, 66,
            139, 1, 142, 79, 5, 132, 2, 174, 227, 106, 143,
            160, 6, 11, 237, 152, 127, 212, 211, 31, 235, 52,
            44, 81, 234, 200, 72, 171, 242, 42, 104, 162, 253,
            58, 206, 204, 181, 112, 14, 86, 8, 12, 118, 18,
            191, 114, 19, 71, 156, 183, 93, 135, 21, 161, 150,
            41, 16, 123, 154, 199, 243, 145, 120, 111, 157, 158,
            178, 177, 50, 117, 25, 61, 255, 53, 138, 126, 109,
            84, 198, 128, 195, 189, 13, 87, 223, 245, 36, 169,
            62, 168, 67, 201, 215, 121, 214, 246, 124, 34, 185,
             3, 224, 15, 236, 222, 122, 148, 176, 188, 220, 232,
            40, 80, 78, 51, 10, 74, 167, 151, 96, 115, 30,
             0, 98, 68, 26, 184, 56, 130, 100, 159, 38, 65,
            173, 69, 70, 146, 39, 94, 85, 47, 140, 163, 165,
            125, 105, 213, 149, 59, 7, 88, 179, 64, 134, 172,
            29, 247, 48, 55, 107, 228, 136, 217, 231, 137, 225,
            27, 131, 73, 76, 63, 248, 254, 141, 83, 170, 144,
            202, 216, 133, 97, 32, 113, 103, 164, 45, 43, 9,
            91, 203, 155, 37, 208, 190, 229, 108, 82, 89, 166,
            116, 210, 230, 244, 180, 192, 209, 102, 175, 194, 57,
            75, 99, 182
           ]

    subInverse = [
            165, 45, 50, 143, 14, 48, 56, 192, 84, 230, 158, 57, 85, 126, 82, 145,
            100, 3, 87, 90, 28, 96, 7, 24, 33, 114, 168, 209, 41, 198, 164, 63, 224,
            39, 141, 12, 130, 234, 174, 180, 154, 99, 73, 229, 66, 228, 21, 183,
            200, 6, 112, 157, 65, 117, 25, 201, 170, 252, 77, 191, 42, 115, 132,
            213, 195, 175, 43, 134, 167, 177, 178, 91, 70, 211, 159, 253, 212, 15,
            156, 47, 155, 67, 239, 217, 121, 182, 83, 127, 193, 240, 35, 231, 37,
            94, 181, 30, 162, 223, 166, 254, 172, 34, 249, 226, 74, 188, 53, 202,
            238, 120, 5, 107, 81, 225, 89, 163, 242, 113, 86, 17, 106, 137, 148,
            101, 140, 187, 119, 60, 123, 40, 171, 210, 49, 222, 196, 95, 204, 207,
            118, 44, 184, 216, 46, 54, 219, 105, 179, 20, 149, 190, 98, 161, 59, 22,
            102, 233, 92, 108, 109, 173, 55, 97, 75, 185, 227, 186, 241, 160, 133,
            131, 218, 71, 197, 176, 51, 250, 150, 111, 110, 194, 246, 80, 255, 93,
            169, 142, 23, 27, 151, 125, 236, 88, 247, 31, 251, 124, 9, 13, 122, 103,
            69, 135, 220, 232, 79, 29, 78, 4, 235, 248, 243, 62, 61, 189, 138, 136,
            221, 205, 11, 19, 152, 2, 147, 128, 144, 208, 36, 52, 203, 237, 244,
            206, 153, 16, 68, 64, 146, 58, 1, 38, 18, 26, 72, 104, 245, 129, 139,
            199, 214, 32, 10, 8, 0, 76, 215, 116
            ]

    #Definition of the AES S-Box for comparison
    subAES = [
            99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171,
            118, 202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156, 164,
            114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113,
            216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39,
            178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227,
            47, 132, 83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74, 76,
            88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80, 60,
            159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16,
            255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61,
            100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20,
            222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98,
            145, 149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244,
            234, 101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198, 232, 221,
            116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14, 97, 53,
            87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155,
            30, 135, 233, 206, 85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104,
            65, 153, 45, 15, 176, 84, 187, 22
            ]

    subAESInverse = [
            82, 9, 106, 213, 48, 54, 165, 56, 191, 64, 163, 158, 129, 243, 215,
            251, 124, 227, 57, 130, 155, 47, 255, 135, 52, 142, 67, 68, 196, 222,
            233, 203, 84, 123, 148, 50, 166, 194, 35, 61, 238, 76, 149, 11, 66, 250,
            195, 78, 8, 46, 161, 102, 40, 217, 36, 178, 118, 91, 162, 73, 109, 139,
            209, 37, 114, 248, 246, 100, 134, 104, 152, 22, 212, 164, 92, 204, 93,
            101, 182, 146, 108, 112, 72, 80, 253, 237, 185, 218, 94, 21, 70, 87,
            167, 141, 157, 132, 144, 216, 171, 0, 140, 188, 211, 10, 247, 228, 88,
            5, 184, 179, 69, 6, 208, 44, 30, 143, 202, 63, 15, 2, 193, 175, 189, 3,
            1, 19, 138, 107, 58, 145, 17, 65, 79, 103, 220, 234, 151, 242, 207, 206,
            240, 180, 230, 115, 150, 172, 116, 34, 231, 173, 53, 133, 226, 249, 55,
            232, 28, 117, 223, 110, 71, 241, 26, 113, 29, 41, 197, 137, 111, 183,
            98, 14, 170, 24, 190, 27, 252, 86, 62, 75, 198, 210, 121, 32, 154, 219,
            192, 254, 120, 205, 90, 244, 31, 221, 168, 51, 136, 7, 199, 49, 177, 18,
            16, 89, 39, 128, 236, 95, 96, 81, 127, 169, 25, 181, 74, 13, 45, 229,
            122, 159, 147, 201, 156, 239, 160, 224, 59, 77, 174, 42, 245, 176, 200,
            235, 187, 60, 131, 83, 153, 97, 23, 43, 4, 126, 186, 119, 214, 38, 225,
            105, 20, 99, 85, 33, 12, 125
            ]

    def __init__(self):
        self.values = [[0] * 8 for _ in range(8)]
  
    def printState(self):
        """Print the state values to the console."""
        for i in self.values:
            print i

    def printStateHex(self):
        """Print the state values as hexadecimal values to the console."""
        for i in range(8):
            out = "%x" % self.getRowHex(i)
            print out.zfill(16)

    def printStatesHex(self, listOfStates):
        """Print a list of states to the console."""
        for row in range(8):
            for state in listOfStates:
                out = "%x" % state.getRowHex(row)
                print out.zfill(16) + " ",
            print ""

    def NumberOfActiveBytes(self):
        """Return the number of non-zero bytes in the state."""
        sumActiveBytes = 0
        for y in range(8):
            for x in range(8):
                if(self.values[x][y] != 0):
                    sumActiveBytes += 1
        return sumActiveBytes 

    def chunks(self, s, n):
        """Produce `n`-character chunks from `s`."""
        chunkList = []
        for start in range(0, len(s), n):
            chunkList.append(s[start:start + n])
        return chunkList
        
    def setStateHex(self, value):
        """Sets the values of the state using hexvalues."""
        # split value in 8 parts
        chunkList = self.chunks(hex(value)[2:-1].zfill(128), 16)
        for row in range(8):
            self.setRowHex(int(chunkList[row], 16), row)
                
    def setRandomNonZeroValue(self, x, y):
        """Set the byte at position (x, y) to a random values in range [1, 255]."""
        self.values[y][x] = random.randint(1, 255)

    def setRandomRow(self, row):
        """Sets a row to random values in range [0, 255]""" 
        self.values[row] = [random.randint(0, 255) for _ in range(8)]
            
    def setRandomValue(self, x, y):
        """Set the byte at position (x, y) to a random values in range [0, 255]."""
        self.values[y][x] = random.randint(0, 255)

    def setValue(self, x, y, value):
        """Set the byte at position (x, y) to value."""
        self.values[y][x] = value

    def getValue(self, x, y):
        """Return the value of the byte at position (x, y)."""
        return self.values[y][x]

    def setFullActive(self):
        """Set all values in the state to 1."""
        self.values = [[1 for _ in range(8)] for _ in range(8)]

    def setActiveRow(self, x):
        """Set all values in a row to 1."""
        self.values[x] = [1 for _ in range(8)]

    def setActiveColumn(self, x):
        """Set all values in a column to 1."""
        for rows in self.values:
            rows[x] = 1

    def Copy(self):
        """Returns a copy of the state."""
        result = GOST()
        for y in range(8):
            for x in range(8):
                result.values[x][y] = self.values[x][y]
        return result

    def AK(self, key):
        """Return the state XORed bytewise with the key."""
        result = GOST()
        for y in range(8):
            for x in range(8):
                result.values[x][y] = self.values[x][y] ^ key.values[x][y]
        return result

    def P(self):
        """Return the state transposed."""
        result = GOST()
        for y in range(8):
            for x in range(8):
                result.values[x][y] = self.values[y][x]
        return result

    def Pinverse(self):
        """Return the state transposed."""
        return self.P()

    def S(self):
        """Return the state after applying the S-Box bytewise."""
        result = GOST()
        result.values = [[self.sub[self.values[j][i]] for i in range(8)] for j in range(8)]
        return result

    def Sinverse(self):
        """Return the state after applying the inverse S-Box bytewise."""
        result = GOST()
        result.values = [[self.subInverse[self.values[j][i]] for i in range(8)] for j in range(8)]
        return result

    def getRow(self, index):
        """Return the row at index."""
        row = []
        for j in range(8):
            row.append(self.values[index][j])
        return row

    def setRow(self, row, index):
        """Set the value of the row at index."""
        for j in range(8):
            self.values[index][j] = row[j]
        return row

    def getRowHex(self, index):
        """Return the row as a 64-bit hex value."""
        row = 0
        for j in range(8):
            row ^= (self.values[index][j]) << ((7 - j) * 8)
        return row

    def setRowHex(self, row, index):
        """Set the row from a 64-bit hex value."""
        for j in range(8):
            self.values[index][j] = (row >> ((7 - j) * 8) & 0xFF)

    def AddConstant(self, index):
        """Return the state after adding the round constant at index."""
        result = GOST()
        if(index >= 0):
            constantState = GOST()
            constantState.setStateHex(self.C[index])
            result = self.AK(constantState)
        return result
        
    def L(self):
        """Return the state after multiplying each row with an 8x8 MDS matrix."""
        result = GOST()
        for rowIndex in range(8):
            # store row in integer
            row = self.getRowHex(rowIndex)
            # perform multiplication
            tmpRow = 0
            for i in range(64):
                if(row & 0x1):
                    tmpRow ^= self.A[i]
                row = row >> 1
            result.setRowHex(tmpRow, rowIndex)
        return result

    def Linverse(self):
        """Return the state after multiplying each row with the inverse of the 8x8 MDS matrix."""
        result = GOST()
        for rowIndex in range(8):
            # store row in integer
            row = self.getRowHex(rowIndex)
            # perform multiplication
            tmpRow = 0
            for i in range(64):
                if(row & 0x1):
                    tmpRow ^= self.Ainverse[i]
                row = row >> 1
            result.setRowHex(tmpRow, rowIndex)
        return result


    def SPL(self):
        """Return the state after applying S, P and L."""
        return self.S().P().L()

    def SPLinverse(self):
        """Return the state after applying Linverse, P and Sinverse."""
        return self.Linverse().P().Sinverse()
