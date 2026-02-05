from dataclasses import dataclass


@dataclass
class NameOptions:
    prefer_alphabetical_tiebreak: bool = True
    return_nd_on_fail: bool = True
    fused_numbering_scheme: str = "iupac2004"
