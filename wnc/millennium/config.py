from cosmap.config.models.sky import Quantity
from cosmap.config.analysis import CosmapAnalysisParameters
from pydantic import Field, BaseModel
from typing import  List

class MsSamplingParameters(BaseModel):
    radius: Quantity
    inner_radius: Quantity
    spacing: int = Field(5, ge=1)
    redshifts: dict[str, float]
    current_field: tuple = (-1,-1)
    fields_to_sample: list[str] = None
    class Config:
        arbitrary_types_allowed = True

class Filters(BaseModel):
    limiting_magnitude: float | list[float]
    filters: dict = {}
    class Config:
        arbitrary_types_allowed = True

class Main(CosmapAnalysisParameters):
    ms_parameters: MsSamplingParameters
    filters: Filters
    weights: List[str]
