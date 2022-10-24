from enum import Enum, unique


@unique
class DomainSuperClass(Enum):
    """
    An Enum representing the superclass of domains supported by RAIChU
    """
    RECOGNITION = 1
    SYNTHESIS = 2
    CARRIER = 3
    TAILORING = 4
    TERMINATION = 5

    @staticmethod
    def from_string(label: str) -> "DomainSuperClass":
        for value in DomainSuperClass:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown domain superclass: {label}")


@unique
class RecognitionDomainType(Enum):
    """
    An Enum representing the types of recognition domains supported by RAIChU
    """
    A = 1
    AT = 2
    DUMMY_A = 3
    DUMMY_AT = 4

    @staticmethod
    def from_string(label: str) -> "RecognitionDomainType":
        for value in RecognitionDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown recognition domain type: {label}")


@unique
class TerminationDomainType(Enum):
    """
    An Enum representing the types of termination domains supported by RAIChU
    """
    TE = 1
    TD = 2
    DUMMY_TE = 3
    DUMMY_TD = 4

    @staticmethod
    def from_string(label: str) -> "TerminationDomainType":
        for value in TerminationDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown termination domain type: {label}")


@unique
class SynthesisDomainType(Enum):
    """
    An Enum representing the types of synthesis domains supported by RAIChU
    """
    C = 1
    KS = 2
    DUMMY_C = 3
    DUMMY_KS = 4

    @staticmethod
    def from_string(label: str) -> "SynthesisDomainType":
        for value in SynthesisDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown synthesis domain type: {label}")


@unique
class CarrierDomainType(Enum):
    """
    An Enum representing the types of carrier domains supported by RAIChU
    """
    PCP = 1
    ACP = 2
    DUMMY_PCP = 3
    DUMMY_ACP = 4

    @staticmethod
    def from_string(label: str) -> "CarrierDomainType":
        for value in CarrierDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown carrier domain type: {label}")


@unique
class TailoringDomainType(Enum):
    """
    An Enum representing the types of tailoring domains supported by RAIChU
    """
    ER = 1
    KR = 2
    DH = 3
    E = 5
    nMT = 6
    UNKNOWN = 7

    @staticmethod
    def from_string(label: str) -> "TailoringDomainType":
        for value in TailoringDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring domain type: {label}")


@unique
class KRDomainSubtype(Enum):
    """
    An Enum representing the subtypes of KR domain supported by RAIChU
    """
    A1 = 1
    A2 = 2
    B1 = 3
    B2 = 4
    C1 = 5
    C2 = 6
    UNKNOWN = 7

    @staticmethod
    def from_string(label: str) -> "KRDomainSubtype":
        for value in KRDomainSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown KR domain subtype: {label}")


@unique
class KSDomainSubtype(Enum):
    """
    An Enum representing the subtypes of KS domain supported by RAIChU
    """
    CIS = 1
    UNKNOWN = 2

    @staticmethod
    def from_string(label: str) -> "KSDomainSubtype":
        for value in KSDomainSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown KS domain subtype: {label}")


@unique
class ERDomainSubtype(Enum):
    """
    An Enum representing the subtypes of KS domain supported by RAIChU
    """
    S = 1
    R = 2
    UNKNOWN = 3

    @staticmethod
    def from_string(label: str) -> "ERDomainSubtype":
        for value in ERDomainSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown ER domain subtype: {label}")
