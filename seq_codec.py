encodePattern = {}
encodePattern[('A', 0)] = 'G'
encodePattern[('A', 1)] = 'C'
encodePattern[('A', 2)] = 'T'

encodePattern[('G', 0)] = 'C'
encodePattern[('G', 1)] = 'T'
encodePattern[('G', 2)] = 'A'

encodePattern[('C', 0)] = 'T'
encodePattern[('C', 1)] = 'A'
encodePattern[('C', 2)] = 'G'

encodePattern[('T', 0)] = 'A'
encodePattern[('T', 1)] = 'G'
encodePattern[('T', 2)] = 'C'


decodePattern = {}
decodePattern[('A', 'A')] = '0'
decodePattern[('A', 'G')] = '0'
decodePattern[('A', 'C')] = '1'
decodePattern[('A', 'T')] = '2'

decodePattern[('G', 'G')] = '0'
decodePattern[('G', 'C')] = '0'
decodePattern[('G', 'T')] = '1'
decodePattern[('G', 'A')] = '2'

decodePattern[('C', 'C')] = '0'
decodePattern[('C', 'T')] = '0'
decodePattern[('C', 'A')] = '1'
decodePattern[('C', 'G')] = '2'

decodePattern[('T', 'T')] = '0'
decodePattern[('T', 'A')] = '0'
decodePattern[('T', 'G')] = '1'
decodePattern[('T', 'C')] = '2'